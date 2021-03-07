#! /bin/bash

celery -A tbprofiler_web.worker worker --loglevel=INFO --concurrency=1
