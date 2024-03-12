#! /usr/bin/python3
##############
# Script for organizing a Cromwell .json input file
# Noah Burget
# 3/12/24
##############
import argparse
def infer_tasks(f, univ_params=[]):
    tasks = []
    with open(f, 'r') as infile:
        for line in infile:
            if '{' in line or '}' in line:
                continue
            key = line.split(':')[0].strip('"').split('.')
            if len(key) < 3:
                univ_params.append(line)
            else:
                if key[1] not in tasks:
                    tasks.append(key[1])
        return tasks, univ_params
def write_organized_file(outf, task_dict, univ_params):
    with open(outf, 'w') as out:
        out.write('{\n')
        print('Writing univrsal parameters')
        for i in univ_params: 
            out.write(i)
        #out.write('\n')
        for tidx, i in enumerate(task_dict):
            print('Writing params for task {}'.format(i))
            for tdidx, j in enumerate(task_dict[i]):
                if tidx == len(task_dict)-1 and tdidx == len(task_dict[i])-1:
                    out.write("{}\n".format(j[:-2]))
                else:
                    out.write(j)
            out.write('\n')
        out.write('}')
        print('Done!')
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, help='Path to input Cromwell JSON')
    args = parser.parse_args()
    tasks, univ_params = infer_tasks(args.input)
    task_dict = {}
    for i in tasks:
        task_dict[i] = []
    #print(task_dict)
    with open(args.input) as f:
        for line in f:
            if '{' in line or '}' in line:
                continue
            key = line.split(':')[0].strip("\t").strip('"').split('.')
            #print(key)
            if len(key) < 3:
                continue
            else:
                if ',' in line[len(line)-2]:
                    task_dict[key[1]].append(line)
                else:
                    task_dict[key[1]].append("{},\n".format(line.split('\n')[0]))
    outf_name = "{}_{}".format(args.input.split('.json')[0],'organized.json')
    write_organized_file(outf_name, task_dict, univ_params)
