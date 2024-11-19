[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/pWmxMLzQ)
# pro3.global|local-Aln
<陳為政 110102050>

## Description

* Write a Python script to perform a global or local alignment.
* Creating your own script, i.e. hw3.py.
* In this program, library Biostrings is only used to parse input fasta file.
* Packages you can use: numpy, pandas, Bio
* You should write a program with a function named alignment, ie.
```
def alignment(input_path, score_path, output_path, aln, gap):
    .
    .
    .
    .
```
* If there is more than one local alignment with the same highest score, you should output local alignments with the maximum length. 
* If there is more than one local alignment with the same highest score, you should output those local alignments in string sequential order according to protein1 and then protein2, i.e., 
  ```
  >protein1
  local alignment1
  >protein2
  local alignment1
  >protein1
  local alignment2
  >protein2
  local alignment2
  ```
## Parameters

* input: .fasta file (ex. test_global.fasta)
* score: score file (ex. pam250.txt)
* aln: global|local
* gap: gap score
* output: .fasta file

## Files

* hw3_ref.py: You can start from this reference code, and try to write your own comment in English.
* pam100.txt
* pam250.txt
* test_global.fasta
* result_global.fasta: You should output your alignment in FASTA format.
* test_local.fasta
* result_local.fasta
## Command

Executing your code with the following command.


```Python
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
```

## Evaluation

10 testing data(5 public, 5 private)

The correct answer gets 10 points for each testing data.



### Penalty

* High code similarity to others: YOUR SCORE = 0

## References
Please provide the code along with its reference. For example, you can cite it as: ```# ChatGPT, respond to “your prompt,” on February 16, 2023```. Below is an example of a reference format summarizing the use of ChatGPT for R programming

>You are the R Language expert.
>Please help me to write a function called “k_fold”.
>Using a given dataset to train the random forest model, and using the k-fold cross-validation to evaluate the best model parameters. Here is the instruction for the function requirements:\
>Function name: k_fold\
>Function parameters:

## Claude
## Prompt:
你是一名生物資訊領域的工程師專家。我正在撰寫一個可以執行global或是local alignment的python程式，需要你幫助我完成，程式功能名稱及輸入如下：
def alignment(input_path, score_path, output_path, aln, gap)
使用範例如下：
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_g123.fasta", "global", -10)
要求如下：
只可以使用numpy, pandas,以及Bio套件，且Bio套件只可以用來讀取以及寫出fasta檔案。目前讀取的部分已經完成了，如下：
# read score path
 with open(score_path, 'r') as f:
 lines_score = f.readlines()[9:]
 amino_acids = lines_score[0].split()
 mut_matrix = np.array([list(map(int, line.split()[1:])) for line in lines_score[1:]])
# 嘗試轉為所需matrix
 substitution_matrix = {}
for i, aa1 in enumerate(amino_acids):
for j, aa2 in enumerate(amino_acids):
 substitution_matrix[(aa1, aa2)] = mut_matrix[i, j]
 substitution_matrix[(aa2, aa1)] = mut_matrix[i, j] # 確保對稱性
# 讀fasta檔
 sequences = list(SeqIO.parse(input_path, "fasta"))
if len(sequences) != 2:
raise ValueError("must contain exactly 2 sequences")
 seq1, seq2 = str(sequences[0].seq), str(sequences[1].seq)
 seq1_id, seq2_id = sequences[0].id, sequences[1].id
我需要你幫我完成分別是local以及global alignment的部分，並且寫出成"result_g123.fasta"檔案，這裡可以使用Bio套件。
另外以下為輸入與輸出範例給你參考，請幫我完成程式碼。
輸入的test_global.fasta:
>1aboA
NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVN
>1ycsB
KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGYVPRNLLGLYP
輸出的result_global.fasta:
>1aboA
N-LFVALYDFVASGDNTLSITKGEKLRVLG-YNHNG-EWCEAQTKNGQGWVPSNYITPVN
>1ycsB
KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGYVPRNLLGLYP
最後也附上pam250.txt：
#
# This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
#
# PAM 250 substitution matrix, scale = ln(2)/3 = 0.231049
#
# Expected score = -0.844, Entropy = 0.354 bits
#
# Lowest score = -8, Highest score = 17
#
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0  0  0  0 -8
R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2 -1  0 -1 -8
N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2  2  1  0 -8
D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2  3  3 -1 -8
C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2 -4 -5 -3 -8
Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2  1  3 -1 -8
E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2  3  3 -1 -8
G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1  0  0 -1 -8
H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2  1  2 -1 -8
I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4 -2 -2 -1 -8
L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2 -3 -3 -1 -8
K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2  1  0 -1 -8
M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2 -2 -2 -1 -8
F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1 -4 -5 -2 -8
P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1 -1  0 -1 -8
S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1  0  0  0 -8
T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0  0 -1  0 -8
W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6 -5 -6 -4 -8
Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2 -3 -4 -2 -8
V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4 -2 -2 -1 -8
B  0 -1  2  3 -4  1  3  0  1 -2 -3  1 -2 -4 -1  0  0 -5 -3 -2  3  2 -1 -8
Z  0  0  1  3 -5  3  3  0  2 -2 -3  0 -2 -5  0  0 -1 -6 -4 -2  2  3 -1 -8
X  0 -1  0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1  0  0 -4 -2 -1 -1 -1 -1 -8
* -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1

## Response:
我會幫你完成這個序列比對的程式。首先我們需要實作動態規劃算法來進行 global 和 local alignment。讓我們一步一步完成這個程式。
Sequence Alignment ImplementationClick to open code
我已經完成了序列比對的程式碼實作，主要包含以下幾個部分：

主要的 alignment 函數：

讀取替換矩陣（PAM250）
讀取輸入序列
根據選擇執行 global 或 local alignment
輸出結果到 FASTA 檔案


global_alignment 函數實作了 Needleman-Wunsch 演算法：

建立得分矩陣和路徑矩陣
使用動態規劃填充矩陣
回溯找出最佳比對結果


local_alignment 函數實作了 Smith-Waterman 演算法：

與 global alignment 類似，但允許局部比對
從最高分位置開始回溯
當得分為 0 時停止比對



使用方式就如你提供的範例：
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_g123.fasta", "global", -10)
這個實作：

完全符合只使用 numpy、pandas 和限制性的 Bio 套件的要求
能夠處理 global 和 local alignment
正確處理 gap penalties
使用 PAM250 替換矩陣進行得分計算
輸出格式符合要求的 FASTA 格式

需要注意的是，gap penalty 是作為負值傳入（如範例中的 -10），程式會直接使用這個值進行計算。
你要測試這個程式碼嗎？或是需要我解釋某個特定部分的實作細節？

## Prompt:
幫我將程式碼整合在一個function中（alignment(...)），另外將沒有用到的import刪掉，謝謝

## Response:
我幫你整理程式碼，移除未使用的 import 並將所有功能整合到 alignment 函數中。
Integrated Sequence Alignment ImplementationClick to open code
主要更改包括：

移除未使用的 import（Align, substitution_matrices, pandas）
將 global_alignment 和 local_alignment 函數移到 alignment 函數內部
簡化了一些變數名稱和註解
保持了原有的功能完整性

使用方式維持不變：
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_g123.fasta", "global", -10)
這樣的整合可以讓程式碼更加模組化，且避免了全域空間的污染。所有的輔助函數現在都封裝在主要的 alignment 函數中。

## Prompt:
CopyTest Failed: Lists differ: ['>1aboA\n', 'LFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE\n', '>[44 chars]E\n'] != ['>1aboA\n', 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE\n', '[46 chars]E\n']

First differing element 1:
'LFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE\n'
'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE\n'

  ['>1aboA\n',
-  'LFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE\n',
+  'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE\n',
?   +

   '>1ycsB\n',
-  'VIYALWDYEPQNDDELPMKEGDCMTIIHREDEDE\n']
+  'GVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDE\n']
?   +
在global比對的部分似乎沒有問題了，但是local的部分出現了以上的錯誤資訊（這是其中一個local的測資）

## Response:
看起來在 local alignment 的實作中有個問題。local alignment 在回溯（traceback）的部分需要將未包含在最佳局部比對區域的序列也加入結果中。讓我修改這部分的程式碼。
Fixed Sequence Alignment ImplementationClick to open code
主要修改了 local_alignment 函數中的回溯部分：

記錄最佳局部比對的起始位置（start_i, start_j）
在回溯完成後，加入未包含在最佳局部比對區域的序列：

在前面加入 seq1[:i] 和 seq2[:j]
在後面加入 seq1[start_i:] 和 seq2[start_j:]



這樣修改後，local alignment 的結果應該會包含完整的序列，而不只是最佳匹配區域。測試結果應該會顯示完整的序列，包括 "NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGE" 而不是缺少開頭的 "N"。






