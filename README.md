# Bioinformatics Project

This project contains multiple tasks (TP2, TP3, TP4, TP5, TP6) related to bioinformatics analysis. Each task is organized in its own folder with relevant data, source code, and documentation.

## Folder Structure

```
bioinfo-main/
    TP2/
        data/
            seq_amont.fa
            sequence.fasta
            sequence.gb
        src/
            utils.py
        test.py
        TP2.md
    TP3/
        data/
            pfm.txt
            seq_amont.fa
        main.py
        scan_pwm.py
        src/
            pwn.py
            utils.py
        TP3.md
    TP4/
        data/
            NM_001100_1000.fa
            NM_001267550_1000.fa
            NM_002469_1000.fa
            NM_002470_1000.fa
            NM_003279_1000.fa
            ...
        main.py
        src/
    TP5/
        data/
        putativeTFBS.py
        README.md
        src/
        test/
    TP6/
        export.json
        main.js
        style.css
        test.json
        viewPutativeTFBS.html
```

## Task Descriptions

### TP2
- **Data Files**: Contains sequence files in various formats (`seq_amont.fa`, `sequence.fasta`, `sequence.gb`).
- **Source Code**: Utility functions in `src/utils.py`.
- **Tests**: `test.py` contains tests for the utility functions.
- **Documentation**: `TP2.md` provides detailed information about the task.

### TP3
- **Data Files**: Contains a position frequency matrix (`pfm.txt`) and a sequence file (`seq_amont.fa`).
- **Source Code**: 
  - `main.py`: Main script for the task.
  - `scan_pwm.py`: Script for scanning sequences with position weight matrices.
  - `src/pwn.py`: Functions for manipulating matrices and scanning sequences.
  - `src/utils.py`: Utility functions.
- **Documentation**: `TP3.md` provides detailed information about the task.

### TP4
- **Data Files**: Contains multiple sequence files (`NM_001100_1000.fa`, `NM_001267550_1000.fa`, etc.).
- **Source Code**: 
  - `main.py`: Main script for the task.
  - `src/`: Directory for additional source code.
- **Documentation**: No specific documentation file provided.

### TP5
- **Data Files**: Directory for data files.
- **Source Code**: 
  - `putativeTFBS.py`: Main script for identifying putative transcription factor binding sites.
  - `src/`: Directory for additional source code.
- **Tests**: `test/` directory contains tests for the source code.
- **Documentation**: `README.md` provides detailed information about installation, usage, and folder structure.

### TP6
- **Data Files**: Contains JSON files (`export.json`, `test.json`).
- **Source Code**: 
  - `main.js`: Main JavaScript file.
  - `style.css`: CSS file for styling.
  - `viewPutativeTFBS.html`: HTML file for viewing putative transcription factor binding sites.
- **Documentation**: No specific documentation file provided.

## Installation and Usage

### TP5 Example
To use the TP5 script, ensure you have Python 3 and Biopython installed. You can install Biopython using pip:

```sh
pip3 install biopython
```

Run the script with the following command:

```sh
python3 putativeTFBS.py -m ./data/pfm.txt -a MA0114 -t 0 -l 1000 -w 40 -s 0 -p 0.1 -mo 5 -o ./export.json NM_001100 NM_002469 NM_002470 NM_003279 NM_003281
```

This command contains 9 parameters followed by the accession numbers of the mRNAs to analyze.

|  Parameter         | Alias | Description | Required |
| :----------------- | :---: | :---------- | :------: |
| --pfm              | -m    | Relative path to the file containing matrices | X |
| --matrix-id        | -a    | Matrix ID to use | X |
| --threshold        | -t    | Score threshold | X |
| --promoter-length  | -l    | Length of the promoter | X |
| --window-size      | -w    | Length of the window | X |
| --start            | -s    | Start position | X |
| --pseudocount      | -p    | Pseudocount value | X |
| --min-occurrences  | -mo   | Minimum occurrences | X |
| --output           | -o    | Output file path | X |

## Authors

- Oudshoorn Gaëtan
- Villegas Navarro Octavio Salvador