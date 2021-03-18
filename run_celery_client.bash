#!/bin/bash
# Run the celery workers

cd ~/similarity_lab
conda activate flaskenv
celery -A app.celery_client worker

