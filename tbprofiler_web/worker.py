from celery import Celery
import subprocess
import json
import sqlite3
import pathogenprofiler as pp
import tbprofiler as tbp
import os
import sys
from flask import Flask,current_app,url_for

import numpy as np

try:
    sys.base_prefix
except:
    sys.base_prefix = getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix))

def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf


def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=app.config['CELERY_RESULT_BACKEND'],
        broker=app.config['CELERY_BROKER_URL']
    )
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    return celery



flask_app = Flask(__name__)
flask_app.config.update(
    CELERY_BROKER_URL='redis://localhost:6379',
    CELERY_RESULT_BACKEND='redis://localhost:6379',
    NEO4J_URI="neo4j://localhost:7687",
    NEO4J_USER="neo4j",
    NEO4J_PASSWORD="test",
)
celery = make_celery(flask_app)

@celery.task
def tbprofiler(fq1,fq2,uniq_id,storage_dir,platform,result_file_dir):

    class Logger(object):
        def __init__(self):
            self.terminal = sys.stdout
            self.log = open("%s/%s.log" % (result_file_dir,uniq_id), "a",buffering=1)

        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

        def flush(self):
            pass

        def revert(self):
            sys.stdout = self.terminal

    # sys.stdout = Logger()
    sys.stderr = Logger()

    conf = get_conf_dict(sys.base_prefix+"/share/tbprofiler/tbdb")
    drug_order = [
        "isoniazid","rifampicin","ethambutol","pyrazinamide","streptomycin",
        "ethionamide","fluoroquinolones","amikacin","capreomycin","kanamycin"
    ]

    if fq1 and fq2:
        fastq_obj = pp.fastq(fq1,fq2)
    elif fq1 and fq2==None:
        fastq_obj = pp.fastq(fq1)
    files_prefix = storage_dir+"/"+uniq_id
    bam_obj = fastq_obj.map_to_ref(
        ref_file=conf["ref"], prefix=files_prefix,sample_name=uniq_id,
        aligner="bwa", platform=platform, threads=4, max_mem="348M"
    )
    bam_file = bam_obj.bam_file

    results = pp.bam_profiler(
        conf=conf, bam_file=bam_file, prefix=files_prefix, platform=platform,
        caller="freebayes", threads=4, no_flagstat=False,
        run_delly = False,samclip=True
    )

    results = tbp.reformat(results, conf, reporting_af=0.1)

    results["id"] = uniq_id
    results["tbprofiler_version"] = tbp._VERSION
    results["pipeline"] = {"mapper":"bwa","variant_caller":"freebayes"}
    outfile = "%s.results.json" % (storage_dir+"/"+uniq_id)





    json_output = "%s/%s.results.json" % (result_file_dir,uniq_id)
    text_output = "%s/%s.results.txt" % (result_file_dir,uniq_id)
    csv_output = "%s/%s.results.csv" % (result_file_dir,uniq_id)
    pdf_output = "%s/%s.results.pdf" % (result_file_dir,uniq_id)
    json.dump(results,open(json_output,"w"))
    tbp.write_pdf(results,conf,pdf_output)
    tbp.write_text(results,conf,text_output)
    tbp.write_csv(results,conf,csv_output)
    pp.run_cmd("samtools view -b -L %s %s > %s/%s.targets.bam" % (conf["bed"], bam_file, result_file_dir, uniq_id))
    pp.run_cmd("samtools index  %s/%s.targets.bam" % (result_file_dir, uniq_id))
    pp.run_cmd("bcftools view %s/%s.targets.csq.vcf.gz > %s/%s.targets.vcf" % (storage_dir,uniq_id,result_file_dir,uniq_id))


    sys.stderr.write("Profiling complete!")
    pp.run_cmd("rm %s/%s*" % (storage_dir,uniq_id))


    sys.stderr.revert()
    return True



@celery.task
def calculate_mutation_stats(gene,variant,pval_cutoff=0.05):
    with flask_app.app_context():
        neo4j_db = get_neo4j_db()

        conf = tbp.get_conf_dict("tbdb")
        lt2drugs = tbp.get_lt2drugs(conf["bed"])
        json_db = json.load(open(conf["json_db"]))


        if gene in json_db and variant in json_db[gene]:
            drugs = list(json_db[gene][variant])
        else:
            drugs = lt2drugs[gene]


        results = []
        for drug in drugs:
            if drug=="para-aminosalicylic_acid": drug = "paraAminosalicylicAcid"
            num_tested = neo4j_db.read("MATCH (s:Sample) WHERE s.%s<>'NA' return s.%s as dst, count(*) as count" % (drug,drug))
            num_tested = {list(x.values())[0]:list(x.values())[1] for x in num_tested }

            variant_data = neo4j_db.read("MATCH (v:Variant {id:'%s_%s'}) <-[:CONTAINS]- (s:Sample) return s.%s as dst, count(*) as count" % (gene,variant,drug))
            variant_data = {list(x.values())[0]:list(x.values())[1] for x in variant_data}
            t = [
                    [0.5, 0.5],
                    [0.5, 0.5]
                ]

            t[0][0] = variant_data.get("1",0.5)
            t[0][1] = variant_data.get("0",0.5)
            t[1][0] = num_tested.get("1",0.5) - variant_data.get("1",0)
            t[1][1] = num_tested.get("0",0.5) - variant_data.get("0",0)

            t2 = sm.stats.Table2x2(np.asarray(t))

            result = {"drug":drug}
            result["odds ratio"] = t2.oddsratio if t[0]!=[0.5,0.5] else "NA"
            result["odds ratio p-value"] = t2.oddsratio_pvalue() if t[0]!=[0.5,0.5] else "NA"
            result["risk ratio"] = t2.riskratio if t[0]!=[0.5,0.5] else "NA"
            result["risk ratio p-value"] = t2.riskratio_pvalue() if t[0]!=[0.5,0.5] else "NA"

            if result["odds ratio"]=="NA":
                result["confidence"] = "indeterminate"
            elif result["odds ratio"]>10 and result["odds ratio p-value"]<pval_cutoff and result["risk ratio"]>1 and result["risk ratio p-value"]<pval_cutoff:
                result["confidence"] = "high"
            elif 5<result["odds ratio"]<=10 and result["odds ratio p-value"]<pval_cutoff and result["risk ratio"]>1 and result["risk ratio p-value"]<pval_cutoff:
                result["confidence"] = "moderate"
            elif 1<result["odds ratio"]<=5 and result["odds ratio p-value"]<pval_cutoff and result["risk ratio"]>1 and result["risk ratio p-value"]<pval_cutoff:
                result["confidence"] = "low"
            elif (result["odds ratio"]<=1 and result["odds ratio p-value"]<pval_cutoff) or (result["risk ratio"]<=1 and result["risk ratio p-value"]<pval_cutoff):
                result["confidence"] = "no_association"
            else:
                result["confidence"] = "indeterminate"

            results.append(result)

        neo4j_db.write(
            "MERGE (v:Variant {id:'%s_%s'})" % (gene,variant),
            "SET v.support = '%s'" % json.dumps(results)
        )
