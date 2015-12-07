#! /bin/python3
# $Id:$
#  --------------------------------------------------------------------------
# dismod_at: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# $begin import_cascade$$ $newlinech #$$
#
# $section Import an IHME Cascade Study$$
#
# $head Syntax$$
# $codei%import_cascade.py %cascade_path% %option_file%$$
#
# $end
# ---------------------------------------------------------------------------
import sys
import os
import csv
import pdb
#
sys.path.append( os.path.join( os.getcwd(), 'python' ) )
import dismod_at
# ---------------------------------------------------------------------------
if sys.argv[0] != 'bin/database2csv.py' :
	msg  = 'bin/database2csv.py: must be executed from its parent directory'
	sys.exit(msg)
#
usage = 'bin/database2csv.py database_file'
if len(sys.argv) != 2 :
	sys.exit(usage)
#
database_file_arg = sys.argv[1]
database_dir      = os.path.split(database_file_arg)[0]
# ----------------------------------------------------------------------------
file_name   = database_file_arg
new         = False
connection  = dismod_at.create_connection(file_name, new)
cmd     = "SELECT * FROM sqlite_master WHERE type='table' AND name='var'"
cursor  = connection.cursor()
result  = cursor.execute(cmd).fetchall()
if len(result) == 0 :
	msg  = 'bin/database2csv.py: must first run init command; i.e.\n'
	msg += '\tdismod_at ' + file_name + ' init'
	sys.exit(msg)
#
cmd = cmd.replace('var', 'sample')
result  = cursor.execute(cmd).fetchall()
have_sample = len(result) > 0
#
cmd = cmd.replace('sample', 'simulate')
result  = cursor.execute(cmd).fetchall()
have_simulate = len(result) > 0
#
cmd = cmd.replace('simulate', 'fit_var')
result  = cursor.execute(cmd).fetchall()
have_fit = len(result) > 0
if have_fit :
	assert len(result) == 1
	cmd = cmd.replace('fit_var', 'fit_data_subset')
	result  = cursor.execute(cmd).fetchall()
	assert len(result) == 1
#
table_data  = dict()
table_list  = [
	'age',
	'covariate',
	'data',
	'data_subset',
	'density',
	'integrand',
	'option',
	'mulcov',
	'node',
	'prior',
	'rate',
	'smooth_grid',
	'time',
	'var',
	'weight'
]
if have_sample :
	table_list.append('sample')
if have_simulate :
	table_list.append('simulate')
if have_fit :
	table_list.append('fit_var')
	table_list.append('fit_data_subset')
for table in table_list :
	table_data[table] = dismod_at.get_table_dict(connection, table)
# ----------------------------------------------------------------------------
# parent_node_id, fit_sample_index
for row in table_data['option'] :
	if row['option_name'] == 'parent_node_id' :
		parent_node_id = int( row['option_value'] )
	#
	if row['option_name'] == 'fit_sample_index' :
		fit_sample_index = row['option_value']
		if fit_sample_index != None :
			sample_index = int(fit_sample_index)
# ----------------------------------------------------------------------------
def convert2output(value_in) :
	if value_in == None :
		value_out = ''
	elif type(value_in) is float :
		value_out = '%13.5g' % value_in
	else :
		value_out = str(value_in)
	return value_out
# ----------------------------------------------------------------------------
def table_lookup(table_name, row_id, column_name) :
	if row_id == None :
		return ''
	value_in = table_data[table_name][row_id][column_name]
	return convert2output(value_in)
# ----------------------------------------------------------------------------
def get_prior_info(row_out, prior_id_dict) :
	extension2name = {'_v':'value_', '_a':'dage_', '_t':'dtime_' }
	for extension in extension2name :
		name      = extension2name[extension]
		key       = name + 'prior_id'
		prior_id  = prior_id_dict[key]
		field_out = 'density' + extension
		if prior_id == None :
			row_out[field_out] = ''
		else :
			density_id   = table_data['prior'][prior_id]['density_id']
			density_name = table_data['density'][density_id]['density_name']
			row_out[field_out] = density_name
		for field_in in [ 'lower', 'upper', 'mean', 'std', 'eta' ] :
			field_out = field_in + extension
			row_out[field_out] = ''
			if prior_id != None :
				value_in = table_data['prior'][prior_id][field_in]
				if field_in == 'eta' :
					if density_name in [ 'uniform', 'gaussian', 'laplace' ] :
						value_in = None
				if field_in == 'std' and density_name == 'uniform' :
						value_in = None
				row_out[field_out] = convert2output( value_in )

# ----------------------------------------------------------------------------
def node_id2child_or_parent(node_id) :
	if node_id == parent_node_id :
		name = table_data['node'][node_id]['node_name']
		return name
	descendant_id = node_id
	while descendant_id != None :
		parent_id  = table_data['node'][descendant_id]['parent']
		if parent_id == parent_node_id :
			name = table_data['node'][descendant_id]['node_name']
			return name
		descendant_id = parent_id
	return None
# ============================================================================
# variable.csv
# ============================================================================
file_name = os.path.join(database_dir, 'variable.csv')
csv_file  = open(file_name, 'w')
#
header = [
	'var_id',
	'var_type',
	'age',
	'time',
	'rate',
	'integrand',
	'covariate',
	'node',
	'fixed',
	'fit_value',
	'fit_res',
	'fit_dage',
	'fit_dtime',
	'sample_value'
]
for extension in ['_v', '_a', '_t' ] :
	for root in ['lower', 'upper', 'mean', 'std', 'eta', 'density' ] :
		header.append( root + extension )
csv_writer = csv.DictWriter(csv_file, fieldnames=header)
csv_writer.writeheader()
var_id  = 0
row_out = dict()
for field in header :
	row_out[field] = ''
for row_in in table_data['var'] :
	row_out['var_id']    = var_id
	row_out['var_type']  = row_in['var_type']
	row_out['age']       = table_lookup('age',  row_in['age_id'], 'age')
	row_out['time']      = table_lookup('time', row_in['time_id'], 'time')
	row_out['rate']      = table_lookup('rate', row_in['rate_id'], 'rate_name')
	row_out['integrand'] = table_lookup(
		'integrand', row_in['integrand_id'], 'integrand_name'
	)
	row_out['covariate'] = table_lookup(
		'covariate', row_in['covariate_id'], 'covariate_name'
	)
	row_out['node'] = table_lookup('node', row_in['node_id'], 'node_name')
	#
	if have_fit :
		row_out['fit_value']= table_lookup('fit_var', var_id, 'variable_value')
		row_out['fit_res']  = table_lookup('fit_var', var_id, 'residual_value')
		row_out['fit_dage'] = table_lookup('fit_var', var_id, 'residual_dage')
		row_out['fit_dtime']= table_lookup('fit_var', var_id, 'residual_dtime')
	if have_sample :
		row_out['sample_value'] = table_lookup('sample', var_id, 'var_value')
	#
	row_out['fixed'] = 'true'
	if row_in['var_type'] == 'rate' :
		if row_in['node_id'] != parent_node_id :
			row_out['fixed'] = 'false'
	#
	smooth_id = row_in['smooth_id']
	if row_in['var_type'].startswith('mulstd_') :
		prior_id_dict = table_data['smooth'][smooth_id]
		get_prior_info(row_out, prior_id_dict)
	else :
		age_id    = row_in['age_id']
		time_id   = row_in['time_id']
		for row in table_data['smooth_grid'] :
			match = row['smooth_id'] == smooth_id
			match = match and row['age_id'] == age_id
			match = match and row['time_id'] == time_id
			if match :
				prior_id_dict = row
				get_prior_info(row_out, prior_id_dict)
	csv_writer.writerow(row_out)
	var_id += 1
# ----------------------------------------------------------------------------
csv_file.close()
# ============================================================================
# data.csv
# ============================================================================
file_name = os.path.join(database_dir, 'data.csv')
csv_file  = open(file_name, 'w')
#
header = [
	'data_id',
	'age_lo',
	'age_up',
	'time_lo',
	'time_up',
	'integrand',
	'weight',
	'hold_out',
	'density',
	'eta',
	'meas_std',
	'meas_value',
	'avgint',
	'residual',
	'node'
]
for row in table_data['covariate'] :
	header.append( row['covariate_name'] )
if fit_sample_index != None :
	index         = header.index('meas_value')
	header[index] = 'sim_value'
csv_writer = csv.DictWriter(csv_file, fieldnames=header)
csv_writer.writeheader()
#
row_out = dict()
for field in header :
	row_out[field] = ''
subset_id  = 0
n_subset   = len( table_data['data_subset'] )
for subset_row in table_data['data_subset'] :
	for field in header :
		row_out[field] = ''
	data_id = subset_row['data_id']
	row_in  = table_data['data'][data_id]
	#
	row_out['data_id']     = data_id
	row_out['age_lo']      = row_in['age_lower']
	row_out['age_up']      = row_in['age_upper']
	row_out['time_lo']     = row_in['time_lower']
	row_out['time_up']     = row_in['time_upper']
	row_out['hold_out']    = row_in['hold_out']
	row_out['meas_std']    = convert2output( row_in['meas_std'] )
	row_out['integrand'] = table_lookup(
		'integrand', row_in['integrand_id'], 'integrand_name'
	)
	row_out['weight'] = table_lookup(
		'weight', row_in['weight_id'], 'weight_name'
	)
	row_out['density'] = table_lookup(
		'density', row_in['density_id'], 'density_name'
	)
	row_out['node'] = node_id2child_or_parent( row_in['node_id'] )
	#
	if row_out['density'] != 'gaussian' :
		row_out['eta']         = table_lookup(
			'integrand', row_in['integrand_id'], 'eta'
		)
	#
	covariate_id = 0
	for row in table_data['covariate'] :
		field_in  = 'x_' + str(covariate_id)
		field_out = row['covariate_name']
		reference = row['reference']
		if row_in[field_in] == None :
			row_out[field_out] = ''
		else :
			row_out[field_out] = convert2output(row_in[field_in] - reference)
		covariate_id += 1
	#
	if have_fit :
		row                 = table_data['fit_data_subset'][subset_id]
		row_out['avgint']   = convert2output( row['avg_integrand'] )
		row_out['residual'] = convert2output( row['weighted_residual'] )
	if fit_sample_index == None :
		row_out['meas_value']  = convert2output( row_in['meas_value'] )
	else :
		if have_simulate :
			simulate_id =  n_subset * sample_index + subset_id
			sim_value = table_data['simulate'][simulate_id]['meas_value']
			row_out['sim_value'] = sim_value
	csv_writer.writerow(row_out)
	subset_id += 1
csv_file.close()
