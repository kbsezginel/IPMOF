# IPMOF Data Analysis Functions
# Date: October 2016
# Author: Kutay B. Sezginel
import os
import math

import yaml
import numpy as np

from ipmof.parameters import read_parameters
from ipmof.interpenetration import regenerate
from tabulate import tabulate


def read_interpenetration_results(results_path):
    """
    Reads interpenetration results from '~/results/S1_S2/results.yaml' file.
    Returns simulation_parameters, structure_info, and summary.
     >>> sim_par, structure_info, summary = read_interpenetration_results(results_path)
    """
    sim_par, structure_info, summary = yaml.load_all(open(results_path, 'r'))
    sim_par = sim_par['simulation_parameters']
    structure_info = structure_info['structure_info']
    summary = summary['summary']
    return sim_par, structure_info, summary


def regenerate_structures(results_dir, export_dir=None, colorify=True, file_format='cif', num=5):
    """
    Regenerates interpenetrated structures from given results.yaml file.
    """
    results_path = os.path.join(results_dir, 'results.yaml')
    with open(results_path, 'r') as r:
        sim_par, structure_info, summary = read_interpenetration_results(results_path)

    if num == 'all':
        num = len(structure_info) - 1
    else:
        num = min(len(structure_info) - 1, num)

    s1_name = structure_info[0]['S1']
    s2_name = structure_info[0]['S2']
    if export_dir is None:
        export_dir = results_dir
    sim_par, sim_dir = read_parameters()

    for s_idx in range(1, num + 1):
        init_coor = structure_info[s_idx]['initial_coordinate']
        rot = structure_info[s_idx]['rotation']
        regenerate(s1_name, s2_name, rot, init_coor, sim_par, sim_dir,
                   export_dir, colorify=colorify, index=s_idx, format=file_format)


def summarize_results(results_dir, summary_path, dir_sep=False, table=True, full=False, format='txt', sortby='structure'):
    """
    Summarize interpenetration results in a given directory and export summary file.
    Optional arguments:
     - dir_sep (True or False):
        must be True if directory separation is used for exporting results

     - table (True or False):
        generate table for min. energy of each MOF combination

     - full (True or False):
        generate table for all combinations (True) or only for combinations
        that resulted in interpenetration (False)

     - sortby ('structure', 'energy', 'collision'):
        sort the table according to number of structures discovered,
        energy of the structure or extended collision
    """
    total_num_of_structures = 0
    table_lines = []
    ip_count = dict(homo=0, hetero=0, homo_total=0, hetero_total=0, no_res=0, no_structure=0, err=0)
    homo_structures = 0
    hetero_structures = 0
    combination_list = []
    error_messages = ''
    no_results = ''

    if os.path.exists(summary_path):
        os.remove(summary_path)
        print('Previous summary removed...')

    if dir_sep:
        for dir_letter in os.listdir(results_dir):
            if os.path.isdir(os.path.join(results_dir, dir_letter)):
                for combination in os.listdir(os.path.join(results_dir, dir_letter)):
                    combination_path = os.path.join(results_dir, dir_letter, combination)
                    combination_list.append(combination_path)
    else:
        combination_list = os.listdir(results_dir)
        combination_list = [os.path.join(results_dir, c) for c in combination_list]

    sim_time_list = []
    for mof_combination in combination_list:
        results_path = os.path.join(mof_combination, 'results.yaml')

        if not os.path.isfile(results_path):
            no_results += '%s\n' % os.path.basename(mof_combination)
            ip_count['no_res'] += 1
        else:
            try:
                sim_par, structure_info, summary = read_interpenetration_results(results_path)
                pass
            except Exception as e:
                ip_count['err'] += 1
                error_messages += '%i. %s -> Error: %s\n' % (ip_count['err'], results_path, e)
                continue
            sim_time = summary['time']
            sim_time_list.append(float(sim_time))
            mof_names = structure_info[0]['S1'] + ' + ' + structure_info[0]['S2']
            num_of_structures = structure_info[0]['Structures']
            total_num_of_structures += num_of_structures

            if structure_info[0]['S1'] == structure_info[0]['S2']:
                ip_count['homo_total'] += 1
            else:
                ip_count['hetero_total'] += 1

            if num_of_structures > 0:
                min_energy = structure_info[1]['energy']
                if structure_info[1]['collision'] is not None:
                    collision = structure_info[1]['collision']['exist']
                else:
                    collision = None
                rotation = structure_info[1]['rotation']
                if structure_info[0]['S1'] == structure_info[0]['S2']:
                    ip_count['homo'] += 1
                    homo_structures += num_of_structures
                else:
                    ip_count['hetero'] += 1
                    hetero_structures += num_of_structures
                if table and not full:
                    table_lines.append([mof_names, num_of_structures, min_energy, collision, rotation, sim_time])
            else:
                ip_count['no_structure'] += 1
                min_energy = collision = rotation = ' - - - '
            if table and full:
                table_lines.append([mof_names, num_of_structures, min_energy, collision, rotation, sim_time])

    ip_count['total'] = ip_count['hetero_total'] + ip_count['homo_total']
    ip_count['structure'] = ip_count['hetero'] + ip_count['homo']
    analysis_text = 'Total MOF combinations:  %i\n' % ip_count['total']
    analysis_text += '\tHetero: %i\tHomo: %i\n' % (ip_count['hetero_total'], ip_count['homo_total'])
    analysis_text += 'Failed jobs:  %i\n' % (ip_count['err'] + ip_count['no_res'])
    analysis_text += '\tNo results file: %i\tCorrupt results file: %i' % (ip_count['no_res'], ip_count['err'])

    if len(sim_time_list) > 0:
        tot_time_s = sum(sim_time_list)
        avg_time = round(tot_time_s / len(sim_time_list), 2)

        minutes, seconds = divmod(tot_time_s, 60)
        hours, minutes = divmod(minutes, 60)
        tot_time = "%d:%02d:%02d" % (hours, minutes, seconds)

        sdr = round(ip_count['structure'] / ip_count['total'] * 100, 2)
        sdr_homo_t = round(ip_count['homo'] / ip_count['total'] * 100, 2)
        sdr_hetero_t = round(ip_count['hetero'] / ip_count['total'] * 100, 2)
        sdpt = round(total_num_of_structures / ip_count['total'], 2)

        if ip_count['homo_total'] == 0:
            sdr_homo_i = '---'
            sdpt_homo_t = '---'
        else:
            sdr_homo_i = str(round(ip_count['homo'] / ip_count['homo_total'] * 100, 2))
            sdpt_homo_t = str(round(homo_structures / ip_count['homo_total'], 2))

        if ip_count['hetero_total'] == 0:
            sdr_hetero_i = '---'
            sdpt_hetero_t = '---'
        else:
            sdr_hetero_i = str(round(ip_count['hetero'] / ip_count['hetero_total'] * 100, 2))
            sdpt_hetero_t = str(round(hetero_structures / ip_count['hetero_total'], 2))

        if ip_count['hetero'] == 0:
            sdpt_hetero_i = '---'
        else:
            sdpt_hetero_i = str(round(hetero_structures / ip_count['hetero'], 2))

        if ip_count['homo'] == 0:
            sdpt_homo_i = '---'
        else:
            sdpt_homo_i = str(round(homo_structures / ip_count['homo'], 2))

        analysis_text += '\nTotal time: %s (%.2f s)\tAverage time: %.2f s\n' % (tot_time, tot_time_s, avg_time)
        analysis_text += '\nMOF combinations with structures:  %i\n' % ip_count['structure']
        analysis_text += '\tHetero:  %i\tHomo: %i\n' % (ip_count['hetero'], ip_count['homo'])
        analysis_text += '\nStructure discovery success (per cent):  %.2f\n\tAmong all trials:\n' % sdr
        analysis_text += '\t\tHetero:  %s\tHomo: %s\n' % (sdr_hetero_t, sdr_homo_t)
        analysis_text += '\tAmong individual trials:\n'
        analysis_text += '\t\tHetero:  %s\tHomo: %s\n' % (sdr_hetero_i, sdr_homo_i)
        analysis_text += '\nTotal number of structures discovered:  %i\n' % total_num_of_structures
        analysis_text += '\tHetero:  %i\tHomo: %i\n' % (hetero_structures, homo_structures)
        analysis_text += '\nStructures discovered per trial (avg):  %.2f\n' % sdpt
        analysis_text += '\tAmong all trials of same type:\n'
        analysis_text += '\t\tHetero:  %s\tHomo: %s\n' % (sdpt_hetero_t, sdpt_homo_t)
        analysis_text += '\tAmong successful trials of same type:\n'
        analysis_text += '\t\tHetero:  %s\tHomo: %s\n\n' % (sdpt_hetero_i, sdpt_homo_i)

    if table and len(table_lines) > 0:
        header = ['MOF Combination', 'Num of Structures', 'Min Energy', 'Collision', 'Rotation', 'Time']
        if sortby == 'structure':
            sort_key = 1
        elif sortby == 'energy':
            sort_key = 2
        elif sortby == 'collision':
            sort_key = 3
        sorted_table = sorted(table_lines, key=lambda x: x[sort_key], reverse=True)

    if format == 'txt' or format == 'both':
        with open(summary_path, 'w') as a:
            a.write(analysis_text)
            if table and len(table_lines) > 0:
                a.write(tabulate(sorted_table, headers=header))
            a.write('\n\nError messages:\n%s\n%s' % ('-' * 100, error_messages))
            a.write('\nResults file(s) could not be found in:\n%s\n%s' % ('-' * 100, no_results))

    if format == 'yaml' or format == 'both' and len(sim_time_list) > 0:
        summary_data = {'tot_combinations': [ip_count['total'], ip_count['hetero_total'], ip_count['homo_total']],
                        'tot_no_res': ip_count['no_res'], 'tot_err': ip_count['err'],
                        'total_time': tot_time_s, 'avg_time': avg_time,
                        'successful_combinations': [ip_count['structure'], ip_count['hetero'], ip_count['homo']],
                        'discovery_rate': [sdr, sdr_hetero_i, sdr_homo_i],
                        'tot_structures': [total_num_of_structures, hetero_structures, homo_structures],
                        'structure_per_trial': [sdpt, sdpt_hetero_t, sdpt_homo_t]
                        }
        if table and len(table_lines) > 0:
            summary_data['table'] = sorted_table

        yaml_summary_path = '%s.yaml' % os.path.splitext(summary_path)[0]
        with open(yaml_summary_path, 'w') as f:
            yaml.dump(summary_data, f)


def get_progress(results_dir, export_dir=None, total=None, dir_sep=False):
    """
    Get progress for number of simulations finished.
    """
    # Choose which directory to export file
    if export_dir is not None and os.path.isdir(export_dir):
        progress_path = os.path.join(export_dir, 'ipmof_progress.txt')
    else:
        progress_path = os.path.join(results_dir, 'ipmof_progress.txt')

    if os.path.exists(progress_path):
        os.remove(progress_path)
        print('Previous progress file removed...')

    table = []
    letter_count = 0
    combination_count = 0
    if dir_sep:
        dir_list = os.listdir(results_dir)
        for dir_letter in dir_list:
            if os.path.isdir(os.path.join(results_dir, dir_letter)):
                letter_count += 1
                combination_list = os.listdir(os.path.join(results_dir, dir_letter))
                combination_count += len(combination_list)
                table.append([dir_letter, len(combination_list), letter_count, combination_count])
    else:
        combination_list = os.listdir(results_dir)
        combination_count = len(combination_list)
        table.append([results_dir, len(combination_list), 1, combination_count])

    header = ['Dir Letter', 'MOF Combinations', 'Dir Count', 'Combination Count']
    sorted_table = sorted(table, key=lambda x: x[1], reverse=True)
    sorted_table = [[i[0], i[1]] for i in sorted_table]

    with open(progress_path, 'w') as p:
        p.write('Total %i combinations in %i directories.\n\n' % (combination_count, letter_count))
        if total is not None:
            p.write('%i combinations left to go.\n\n' % int(total - combination_count))
        p.write(tabulate(table, headers=header))
        p.write('\n\nTable with directories sorted according to the number of combinations:\n\n')
        p.write(tabulate(sorted_table, headers=header[:2]))


def get_scatter_data(results_dir, vf_list, tag='', export_dir=None, error=True):
    """ Get scatter data from a directory of results.yaml files
        kwargs:
         - error = True / False
           Error information can be printed for missing / corrupt results.
    """
    if export_dir is None:
        export_dir = os.getcwd()

    scatter_data = []
    error_messages = ''

    combination_list = os.listdir(results_dir)
    combination_list = [os.path.join(results_dir, i) for i in combination_list]

    for mof_combination in combination_list:
        results_path = os.path.join(mof_combination, 'results.yaml')
        if not os.path.isfile(results_path):
            error_messages += '%s -> No results file\n' % results_path
        else:
            try:
                sim_par, structure_info, summary = read_interpenetration_results(results_path)
                pass
            except Exception as e:
                error_messages += '%s -> Error: %s\n' % (results_path, e)
                continue
            mof1 = structure_info[0]['S1']
            mof2 = structure_info[0]['S2']
            mof1_vf = vf_list['vf'][vf_list['mof'].index(mof1)]
            mof2_vf = vf_list['vf'][vf_list['mof'].index(mof2)]
            num_of_structures = structure_info[0]['Structures']

            scatter_data.append([mof1, mof2, mof1_vf, mof2_vf, num_of_structures])

    scatter_data_path = os.path.join(export_dir, 'scatter_data_%s' % tag)
    scatter_numpy = np.array(scatter_data)
    np.save(scatter_data_path, scatter_numpy)

    if error:
        error_log_path = os.path.join(export_dir, 'scatter_error_%s.txt' % tag)
        with open(error_log_path, 'w') as f:
            f.write(error_messages)
