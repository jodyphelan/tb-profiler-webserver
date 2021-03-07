from flask import (
	Blueprint, flash, g, redirect, render_template, request, url_for, Response
)
from werkzeug.exceptions import abort
import json
# from tbprofiler_web.auth import login_required
import tbprofiler as tbp
from flask import current_app as app
bp = Blueprint('results', __name__)
import os

def get_result(sample_id):
	return json.load(open(app.config["APP_ROOT"]+url_for('static', filename='results/') + sample_id + ".results.json"))

@bp.route('/results/json/<sample_id>',methods=('GET', 'POST'))
def run_result_json(sample_id):
	log_file = app.config["APP_ROOT"]+url_for('static', filename='results/') + sample_id + ".log"
	if os.path.isfile(log_file)==False:
		return {"status":"Invalid_ID","result":None}
	progress = check_progress(log_file)
	if progress!="Completed":
			return {"status":progress,"result":None}

	results = get_result(sample_id)
	return {"status":"OK","result":results}



@bp.route('/results/<sample_id>')
def run_result(sample_id):


	log_file = app.config["APP_ROOT"]+url_for('static', filename='results/') + sample_id + ".log"
	if os.path.isfile(log_file)==False:
		flash("Error! Result with ID:%s doesn't exist" % sample_id)
		return redirect(url_for('home.index'))
	progress = check_progress(log_file)
	print(progress)
	if progress!="Completed":
		log_text = open(log_file).read().replace(app.config["UPLOAD_FOLDER"]+"/","")
		return render_template('results/run_result.html',result = None,sample_id=sample_id,progress = progress,log_text=log_text)

	results = get_result(sample_id)

	for var in results["dr_variants"]:
		print(var["drugs"])
		var["drugs"] = ", ".join([d["drug"] for d in var["drugs"]])

	bam_found = os.path.isfile(app.config["APP_ROOT"]+url_for('static', filename='results/') + sample_id + ".targets.bam")
	return render_template('results/run_result.html',result = results, bam_found = bam_found, sample_id=sample_id)


def check_progress(filename):
	progress = "Uploaded"
	text = open(filename).read()
	if "bwa mem" in text:
		progress = "Mapping"
	if "samtools fixmate"  in text:
		progress = "Bam sorting"
	if "samclip" in text:
		progress = "Variant calling"
	if "bcftools csq" in text:
		progress = "Variant annotation"
	if "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD]" in text:
		progress = "Lineage determination"
	if "Profiling complete!" in text:
		progress = "Completed"
	return progress
