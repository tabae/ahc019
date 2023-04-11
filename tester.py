#coding: utf-8
#requirements: python3.8

import os
import sys
import time
import logging
import datetime
import subprocess

min_seed = 0
max_seed = 99

def compile(src_code, dest_exe):
    cmd = "g++ -std=c++17 -O2 {} -o {}".format(src_code, dest_exe) 
    logging.info("start compile: {}".format(cmd))
    if os.path.isfile(dest_exe):
        os.remove(dest_exe)
    subprocess.run(cmd.split(), check=True)
    if not os.path.isfile(dest_exe):
        logging.error("failed to compile: {}".format(cmd))
        sys.exit(1)
    logging.info("finish compile")

def run(exe_file, log_raw_file):
    logging.info("start run")
    os.chdir('tools')
    logging.info("change directory: {}".format(os.getcwd()))
    for cur_seed in range(min_seed, max_seed+1):
        cur_seed_str = str(cur_seed).zfill(4)
        input_file = os.path.join(os.getcwd(), 'in', '{}.txt'.format(cur_seed_str))
        output_file = os.path.join(os.getcwd(), 'out', '{}.txt'.format(cur_seed_str))
        with open(input_file, mode='r') as f_in, open(output_file, mode='w') as f_out, open(log_raw_file, mode='a') as f_log:
            cmd = "{}".format(exe_file)
            logging.info("run: {}".format(cmd))
            f_log.write('seed = {}\n'.format(cur_seed))
            start_time = time.time()
            subprocess.run(cmd.split(), check=True, stdin=f_in, stdout=f_out, stderr=f_log)
            end_time = time.time()
            f_log.write('elapsed time = {}\n'.format(end_time - start_time))
            cmd = "cargo run --release --bin vis {} {}".format(input_file, output_file)
            logging.info("run: {}".format(cmd))
            subprocess.run(cmd.split(), check=True, stdout=f_log, stderr=f_log)
            f_log.write('==============\n')
    os.chdir('../')
    logging.info("change directory: {}".format(os.getcwd()))
    logging.info("finish run")
        
def summarize(log_raw_file, score_csv_file, log_summary_file, now):
    logging.info('start summarize')
    seeds = []
    myscores = []
    scores = []
    elapsed_times = []
    inputDs = []
    with open(log_raw_file, mode='r') as f_in:
        sc = None; msc = None; elt = None; sd = None; id = None
        for line in f_in:
            if(line.startswith('seed')):
                sd = int(line.split('=')[1])
            if(line.startswith('elapsed time')):
                elt = float(line.split('=')[1])
            if(line.startswith('inputs.D')):
                id = int(line.split('=')[1])
            if(line.startswith('MyScore')):
                msc = int(line.split('=')[1])
            if(line.startswith('Score')):
                sc = int(line.split('=')[1])
            if(line.startswith('====')):
                seeds.append(sd)
                elapsed_times.append(elt)
                myscores.append(msc)
                scores.append(sc)
                inputDs.append(id)
                sc = None; msc = None; elt = None; sd = None; id =None
    if(len(myscores) != len(scores)):
        logging.warning('len(myscores) {} is not same as len(scores) {}'.format(len(myscores), len(scores)))
    else:
        for i in range(len(scores)):
            if(scores[i] != myscores[i]):
                logging.warning('scores[{}] {} != myscores[{}] {}'.format(i, scores[i], i, myscores[i]))
                break
    if(len(scores) != max_seed - min_seed + 1):
        logging.warning('len(scores) {} is not same as num of testcases {}'.format(len(scores), max_seed - min_seed + 1))
    logging.info('longest elapsed time = {:.2f}s'.format(max([0 if ele is None else ele for ele in elapsed_times])))
    logging.info('sum of scores = {:,}'.format(sum([0 if ele is None else ele for ele in scores])))
    cmd = "git show --no-patch --format='%h'"
    res = subprocess.run(cmd.split(), check=True, stdout=subprocess.PIPE)
    githash = res.stdout.decode(sys.getfilesystemencoding()).rstrip().replace('\'', '')
    cmd = "git show --no-patch --format='%s'"
    res = subprocess.run(cmd.split(), check=True, stdout=subprocess.PIPE)
    gitmsg = res.stdout.decode(sys.getfilesystemencoding()).rstrip()
    logging.info('write csv file {}'.format(score_csv_file))
    with open(score_csv_file, mode='w') as f_out:
        f_out.write('#seed, inputsD score, myscore, elapsed time(s)\n')
        for i in range(len(seeds)):
            f_out.write('{}, {}, {}, {}, {}\n'.format(seeds[i], inputDs[i], scores[i], myscores[i], elapsed_times[i]))
    logging.info('write summary log {}'.format(log_summary_file))
    with open(log_summary_file, mode='a') as f_out:
        f_out.write('===================================\n')
        f_out.write('Running date   : {}\n'.format(now.strftime('%y/%m/%d %H:%M:%S')))
        f_out.write('Git hash       : {}\n'.format(githash))
        f_out.write('Git message    : {}\n'.format(gitmsg))
        f_out.write('Sum of score   : {:,}\n'.format(sum([0 if ele is None else ele for ele in scores])))
        f_out.write('Max elaps time : {:.2f}s\n'.format(max([0 if ele is None else ele for ele in elapsed_times])))
        f_out.write('Seeds          : {}-{}\n'.format(min_seed, max_seed))
        f_out.write('Invalid case   : {}/{}\n'.format(scores.count(None), max_seed-min_seed+1))
        f_out.write('Enable         : True\n')
    logging.info('finish summarize')

def main():
    
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.DEBUG,
        format="[%(levelname)s] %(message)s",
    )
    
    now = datetime.datetime.now()
    now_str = now.strftime('%y%m%d-%H:%M:%S')

    # define log path
    log_summary_file = os.path.join(os.getcwd(), 'logs', 'log_summary_file.log')
    log_dir = os.path.join(os.getcwd(), 'logs', now_str)
    log_raw_file = os.path.join(log_dir, 'raw.log')
    score_csv_file = os.path.join(log_dir, 'score.csv')

    # define source file and executive file
    source_file = os.path.join(os.getcwd(), 'main.cpp')
    executive_file = os.path.join(os.getcwd(), 'tools', 'solver')

    # make directory
    if not os.path.isdir(os.path.join(os.getcwd(), 'tools', 'out')):
        os.makedirs(os.path.join(os.getcwd(), 'tools', 'out'))
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logging.info('log file is in {}'.format(log_raw_file))
    
    compile(source_file, executive_file)
    run(executive_file, log_raw_file)
    summarize(log_raw_file, score_csv_file, log_summary_file, now)
    

if __name__ == '__main__':
    main()

