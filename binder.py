import xlrd, os
import shutil
from os import listdir
from os.path import isfile, join

split_threshold = 119.5
column_index = 5
xl_path = r'D:\ASU Classes\Spring17\SML\brain-analysis-for-gender-classification\datasets\set1'
input_data_path = r"D:\ASU Classes\Spring17\SML\brain-analysis-for-gender-classification\datasets\set1\normalized\merged"
split_dir_path_1 = r"D:\ASU Classes\Spring17\SML\brain-analysis-for-gender-classification\datasets\fsiq\low"
split_dir_path_2 = r"D:\ASU Classes\Spring17\SML\brain-analysis-for-gender-classification\datasets\fsiq\high"

def main():
    workbook = xlrd.open_workbook(os.path.join(xl_path, 'metainfo.xls'))
    sheet = workbook.sheet_by_index(0)
    file_names = [f for f in listdir(normalized_data_path) if isfile(join(input_data_path, f))]
    for i in range(1, sheet.nrows):
        for file_name in file_names:
            if sheet.cell(i, 0).value in file_name:
                if sheet.cell(i, column_index).value > split_threshold:
                    shutil.copy(join(input_data_path, file_name), split_dir_path_2)
                else:
                    shutil.copy(join(input_data_path, file_name), split_dir_path_1)
                break

if __name__ == "__main__":
    main()
