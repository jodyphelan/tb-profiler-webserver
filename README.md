# tb-profiler-webserver

This repository hosts the code to deploy a webserver to wrap around the function of [TB-Profiler](https://github.com/jodyphelan/TBProfiler/). Updates will follow soon!

## Installation
Installation requires tb-profiler, flask, celery and redis.
To run it on your local machine:
```
# Install libraries
conda install python=3.7 flask redis celery statsmodels bwa samtools bcftools parallel freebayes gatk4 bedtools samclip delly
pip install neo4j redis tqdm

python setup.py install

# Run flask
export FLASK_APP=tbprofiler_web
export FLASK_ENV=development
flask run

# Run rabbit-mq server
rabbitmq-server

# Run celery
celery -A tbprofiler_web.worker worker --loglevel=info --concurrency=1
```
