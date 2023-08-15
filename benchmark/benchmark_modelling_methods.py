#!/usr/bin/env python3

import datetime
import logging
import os
import subprocess
import sys
import time
import argparse


def parse_chromosomes(input_str):
    if input_str == "autosomal":
        return list(range(1, 22 + 1))

    if '-' in input_str:
        start_chr, end_chr = input_str.split('-')
        return list(range(int(start_chr.replace('chr', '')), int(end_chr.replace('chr', '')) + 1))

    return [int(c.replace('chr', '')) for c in input_str.split(',')]

class Tester:
    def __init__(self, stg_filename, output_dir_name, method='cudaMMC', num_iter=1, num_ensemble=0,
                 which_chromosomes=(1, 14, 21)):
        self.stg_init_file = stg_filename
        self.method = method
        self.num_iter = num_iter
        self.num_ensemble = num_ensemble
        self.which_chromosomes = which_chromosomes
        self.output_dirname = output_dir_name

        if method not in ['cudaMMC', '3dnome']:
            raise ValueError('Wrong input method!')

        self.logs_path = f'logs/{self.method}_logs_tester/{self.output_dirname}_logs'
        self.output_path = f'models/{self.method}_models_tester/{self.output_dirname}_models'

    def __ensure_directories_exist(self):
        os.makedirs(self.logs_path, exist_ok=True)
        os.makedirs(self.output_path, exist_ok=True)

    def __build_command(self, chromosome):
        command = [f'./{self.method}', '-s', self.stg_init_file, '-c', f'chr{chromosome}', '-o',
                   f'{self.output_path}/chr{chromosome}_']

        if self.num_ensemble != 0:
            command.extend(['-m', str(self.num_ensemble)])

        return command

    def test(self):
        self.__ensure_directories_exist()

        for i in range(self.num_iter):
            for c in self.which_chromosomes:
                logging.info(f'Executing {self.method} for chromosome {c}_{i}')
                start = time.perf_counter()

                try:
                    command = self.__build_command(c)
                    logging.info(" ".join(command))
                    subprocess.run(command, check=True, capture_output=False)

                except subprocess.CalledProcessError as e:
                    logging.error(f"{self.method} execution failed! Writing process output to dump.txt")
                    with open('dump.txt', 'w') as dump:
                        dump.write(e.stdout.decode("utf-8"))
                    exit(1)

                duration = time.perf_counter() - start
                result_file_name = f'{self.logs_path}/{c}_{datetime.datetime.today()}.txt'
                logging.info(f"Saving result to the {result_file_name}.")

                with open(result_file_name, 'w') as result_file:
                    result_file.write(str(duration))
if __name__ == "__main__":
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser(description="Run the testing procedure.")

    parser.add_argument("-c", "--config_file", type=str, required=True, help="Configuration file path.")
    parser.add_argument("-o", "--output_dirname", type=str, required=True, help="Output directory name.")
    parser.add_argument("-m", "--method_name", type=str, choices=['cudaMMC', '3dnome'], required=True,
                        help="Method name, either 'cudaMMC' or '3dnome'.")
    parser.add_argument("-i", "--iteration_number", type=int, required=True, help="Number of iterations.")
    parser.add_argument("-e", "--ensemble_number", type=int, default=1, help="Number of ensembles. Defaults to 1.")
    parser.add_argument("-d", "--which_chromosome", type=str, required=True,
                        help="Chromosomes to consider, can be 'autosomal', a range (e.g., 'chr10-chr15'), comma-separated (e.g., 'chr1,chr14,chr21'), or single chromosome (e.g., 'chr21').")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    which_chromosome = parse_chromosomes(args.which_chromosome)

    tester = Tester(
        args.config_file,
        args.output_dirname,
        method=args.method_name,
        num_iter=args.iteration_number,
        num_ensemble=args.ensemble_number,
        which_chromosomes=which_chromosome
    )

    tester.test()