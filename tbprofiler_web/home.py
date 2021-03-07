from flask import (
	Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort
import json
#from aseantb.auth import login_required
from flask import current_app as app

bp = Blueprint('home', __name__)


@bp.route('/')
def index():
	return render_template('home/index.html')

@bp.route('/robots.txt')
def robots():
	# neo4j = get_neo4j_db()
	return open(app.config["APP_ROOT"]+url_for('static',filename='robots.txt')).read().replace("\n","<br>")
