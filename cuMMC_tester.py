#!/usr/bin/env python3

import datetime
import logging
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

class Tester:
    def __init__(self):
        # IAM user credentials that can only read and list objects from roszcdam-gm12878-test s3 bucket
        # and write to roszcdam-gm12878-test/results/, so you can try to mine crypto, but it won't happen :)
        ACCESS_KEY_ID = "AKIAUIO56OEHNKBYU7GH"
        SECRET_ACCESS_KEY = "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"

        self.DATA_FOLDER = 'data_GM12878'
        self.BUCKET_NAME = 'roszcdam-gm12878-test'
        self.REPO_NAME = 'totally_not_cummc'
        self.REPO_URL = f'https://github.com/gprabowski/{self.REPO_NAME}.git'
        self.s3_client = boto3.client(
            's3',
            aws_access_key_id=ACCESS_KEY_ID,
            aws_secret_access_key=SECRET_ACCESS_KEY
        )

    def test(self):
        self.__download_data_from_s3()
        self.__clone_repo()
        self.__compile()
        self.__run()

    def __download_data_from_s3(self):
        logging.info("Attempting to download the data from s3...")

        try:
            os.mkdir(self.DATA_FOLDER)
        except FileExistsError:
            logging.info('Data directory already exists, skipping download...')
            return

        response = self.s3_client.list_objects_v2(Bucket=self.BUCKET_NAME)

        for object in response["Contents"]:
            object_name = object["Key"]

            if 'results' in object_name:
                continue

            logging.info(f'Downloading {object_name}')
            self.s3_client.download_file(
                self.BUCKET_NAME,
                object_name,
                f'{self.DATA_FOLDER}/{object_name}'
            )

    def __clone_repo(self):
        if os.path.isdir(f'./{self.REPO_NAME}'):
            logging.info("Code is already there, skipping.")
            os.chdir(self.REPO_NAME)
            return

        logging.info("Cloning the repo...")

        try:
            subprocess.run(['git', 'clone', self.REPO_URL], check=True)
        except subprocess.CalledProcessError:
            logging.error("Something went wrong while cloning the repo")
            exit(1)

        os.chdir(self.REPO_NAME)

    def __compile(self):
        logging.info("Compiling...")
        try:
            subprocess.run(['make', 'clean'], check=True)
            subprocess.run(['make'], check=True)
        except subprocess.CalledProcessError:
            logging.error("Compilation failed!")
            exit(1)

    def __run(self):
        for c in [1, 14, 21]:
            for setting in ['cache', 'nocache']:
                logging.info(f'Executing cuMMC for chromosome {c} with {setting}')
                start = time.perf_counter()

                try:
                    logging.info(f'./cuMMC -s stg_gpu.ini -c chr{str(c)} -o ./{c}/')
                    subprocess.run(
                        [ './cuMMC', '-s', 'stg_gpu.ini', '-c', 'chr' + str(c), '-o', f'./chr{c}/' ],
                        check=True,
                        capture_output=True
                    )
                except subprocess.CalledProcessError as e:
                    dump_file = 'dump.txt'
                    logging.error(f"cuMMC execution failed! Writing process output to {dump_file}")
                    with open(dump_file, 'w') as dump:
                        dump.write(e.stdout.decode("utf-8"))
                    exit(1)

                end = time.perf_counter()
                result = end - start
                result_file_name = f'{c}_{setting}_{str(datetime.datetime.today())}.txt'

                logging.info(f"Saving result to the {result_file_name}.")
                with open(result_file_name, 'w') as result_file:
                    result_file.write(str(result))

                logging.info(f"Uploading {result_file_name} to s3.")
                try:
                    self.s3_client.upload_file(
                        result_file_name,
                        self.BUCKET_NAME,
                        f'results-{datetime.date.today()}/{result_file_name}'
                    )
                except ClientError:
                    logging.error(f"File {result_file_name} failed to upload to s3!")

if __name__ == "__main__":
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s - %(message)s'
    )
    tester = Tester()
    tester.test()
