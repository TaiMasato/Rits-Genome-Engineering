# Rits-genome-engineering
For-programming-at-Rits-Sugano-Group

本ツールはCRISPR/Casシステムで使用するgRNAを遺伝子の中央および両端に設計する。

深尾研Sugano-Group用のプログラミングを共有するためのGithubリポジトリです。
otiai10さんもアドバイザーとして2017.7.16現在は参加されてます。

# 構成

```
├── NGS_sORF      #NGS解析の時に使ったやつ
│     ├── assets
│     └── src
├── README.md
├── R_script
├── assets        #inputなどが入る
│     ├── index   #bowtieに使うファイルが格納される
│     └── results #結果が格納される
└── src           #ソースコードが格納される

```

# 環境構築 for Mac (本ツールインストール前に)
## 必要ツールおよびライブラリ
>python3   
>numpy 
>bowtie

## インストール
>command line tools for Xcodeのインストール
>```
>	$ xcode-select --install
>```
>Homebrewのインストール
>```
>	$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
>```
>非公式パッケージをインストールするために
>```
>	$ brew tap brewsci/bio
>```
>	または
>```
>	$ brew tap brewsci/science
>```
>Bowtieのインストール
>```
>	$ brew install bowtie
>```
>Bowtieにインデックスを作成
>	assets内に作りたいindexに対応した生物種のfastaをブチ込む。
>```
>	$ bowtie-build -f <fastaファイル名> <命名したいインデックス名>
>```
>	全部で６個のファイルが作られる。
>python3のインストール
>```
>	$ brew install python3
>	$ which python3
>```
>numpyのインストール
>```
> $ brew install numpy
>```
>本ツールのインストール
>```
> $ https://github.com/shigeohu/Rits-genome-engineering.git
>```

# ワークフロー
<img width="1235" alt="スクリーンショット 2019-11-05 14 40 14" src="https://user-images.githubusercontent.com/30248550/68181318-39d93380-ffda-11e9-8d9f-6e656a9c8786.png">
本ツールのアルゴリズムは各生物種のゲノム配列 (Fasta)、遺伝子情報やRNAseqデータ (GFF3)、遺伝子名リスト (txt)、任意のPAM配列およびgRNA設計位置などの入力情報を取得しgRNAの自動設計、評価を一括で行う。マッピングツールであるbowtie (Langmead et al., 2009)を用いて設計されたgRNAのオフターゲットを高速に探索する。自動評価後、各遺伝子の破壊に最も適したオフターゲットが少ないgRNAはGFF3形式のファイルとして出力される。


# Usage
## コマンド
>Framshift: 
>```
>$ python3 main_script.py <input1> <input2> <input_file1> <input_file2> <input_file3> <input3> <input4> <input5> <input6> <input7> <input8>
>```
>Large Deletion: 
>```
>$ python3 main_script.py <input1> <input2> <input_file1> <input_file2> <input_file3> <input3> <input4> <input5> <input8>
>```

## 各種input説明
> ＜input1＞  gRNA designing position: gRNAの設計位置を指定。frameshit位置に設計→ [fs], large deletion位置に設計→ [ld]  
> ＜input2＞  input species number  
>　　　[0 -> *A_thaliana*    ]  
>　　　[1 -> *M_polymorpha*    ]  
>　　　[2 -> *N_Benthamiana*    ]  
> ＜input_file1＞  query input files are fasta format <.fasta/.fa>  
> ＜input_file2＞  query input files are gff3 format <.gff/.gff3>  
> ＜input_file3＞  query input files are gene id list file <.txt>  
> ＜input3＞  PAM: 任意のPAMを入力 ex) NGG  
> ＜input4＞  GC含量を考慮するか否か→ [y/n]。yの場合、GC含量は25~75%保持しているgRNA以外を除外。  
> ＜input5＞  poly_Tを考慮するか否か→ [y/n]。yの場合、Poly_Tは連続するTが4つ以上あるものを除外する。  
> ＜input6＞  真ん中からの取り出しを有効にするか否か→ [y/n]。  
> ＜input7＞  splicing variantを考慮するか否か→ [y/n]。  
> ＜input8＞  Number of gRNAs per gene to extract→ [<int>/n]。最終的に各設計領域から取り出すgRNAの数を指定する。nの場合、各領域に設計されたgRNAの最大数(直前の操作に依存)を取り出す。

# 各ファイルの概要
 0. :integrate_src_control.py  
> 下記のコードを一括で動かす。全ての出力ファイルはここで作製される。全てのinput情報はここで入力する。出力ファイルもここで生成される。
 1. :message.py
> input情報が正しくない場合に表示されるUsageが書かれたコード。
 2. :def_process_input_files.py
> 入力された情報をもとにfasta,gff,gene_idリストを展開する。
 3. :def_1v2_design_all_spacer_frameshift.py
> CDS領域にgRNAを設計する。 fastaとgffを元にCDS領域を探し出し、そこからPAM配列の下流20bpをgRNAとして設計する。設計したgRNAはgff形式ものとfasta形式の2種類が出力される。
 4. :def_1v3_design_all_spacer_large_deletion.py
> 遺伝子の両端にgRNAを設計する。fastaとgffを元に遺伝子間領域を探し出し、そこからPAM配列の下流20bpをgRNAとして設計する。設計したgRNAはgff形式ものとfasta形式の2種類が出力される。
 5. :def_bowtie_process.py
> オフターゲットの探索。fasta形式のgRNAファイルをbowtieに処理させ、オフターゲットファイルを出力させる。
 6. :def_2_pam_check.py
> gRNA類似配列の上流にはPAM配列があるのかを調べる。
 7. :def_3_seed_check.py
> オフターゲット内のミスマッチがSeed領域かnon-Seed領域のどの位置に存在するかを調べる。gRNA上流12bpがSeed領域、残り8bpがSeed領域。
 8. :def_4_offtarget_count.py
> それぞれのgRNAに存在するオフターゲットの数をカウントする。
 9. :def_5_add_OT_num.py
> オフターゲット情報をgff形式gRNAファイルに追加する。
10. :def_6_join_bowtie_summary_2query.py
> フォーマット変換１
11. :def_7_Make_Unique_gRNA_list_with_Annotation.py
> フォーマット変換２
12. :def_8_Restrict_gRNA.py
> (option1) GC含量が25~75%のものおよびPolyTが存在しないgRNAを抽出する。
13. :def_9_extract_gRNA_at_central.py
> (option2) CDSの中央領域からgRNAを10本取り出す。フレームシフト用gRNA設計時のみ使用。
14. :def_10_pick_out_gRNA.py
>  (option3) 各遺伝子ごとに自分が指定した本数のgRNAを取り出す。取り出す順番は最も特異的なgRNA順とする。指定した本数のgRNAが取り出せない場合その遺伝子に存在する最大量のgRNAを取り出す。



