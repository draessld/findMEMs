#   Find All Maximal Exat Matches in the collection of texts
- small project, that build a suffix tree over text collection given in fasta file or as a simple text delimited by special character '$'.
- returns output in the following format:
```
>GTGGGGG        73us    1
GTGGGGG:1       5{0}

>CCCCATATCCCCCC 756us   3
ATCCCCCC:5      74{0}   49{0}   134{0}  174{0}  142{0}
TATCCCC:1       105{0}
CCCCATA:3       89{0}   120{0}  79{0}

>ATCCCCCG       286us   1
ATCCCCCG:1      195{0}
```

##  Compilation
- you need to have sdsl library installed
- please change in Makefile to include your path to sdsl
- run make and have fun

##  Usage

To build index from given text, txt format should be in format of "text1$text2$...$"
```shell
./build test_data/input.txt
./build test_data/input.fasta
./build test_data/input.fa 
```

To locate all MEMs of patterns with respect to the texts.
```shell
./locate test_data/input.txt test_data/patterns.txt
./locate test_data/input.fasta test_data/patterns.fasta
./locate test_data/input.fa test_data/patterns.fa
```

```shell
./locate test_data/input.txt CCCCATATCCCCCCA
./locate test_data/input.fasta CCCCATATCCCCCCA
./locate test_data/input.fa CCCCATATCCCCCCA
```