# tb-profiler-webserver

This repository hosts the code to deploy a webserver to wrap around the function of [TB-Profiler](https://github.com/jodyphelan/TBProfiler/). Updates will follow soon!

## Installation
Installation requires tb-profiler, flask, celery and redis.
To run it on your local machine:
```
# Install libraries
conda install python=3.7 flask redis celery tb-profiler
pip install neo4j redis

git clone https://github.com/jodyphelan/tb-profiler-webserver.git
de tb-profiler-webserver
python setup.py install

# Run flask
bin/flash

# Run redis
redis-server

# Run celery
bin/celery
```
