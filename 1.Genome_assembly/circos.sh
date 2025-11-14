#!/usr/bin/env bash
perl ../scripts/geneid_to_mRNAid.pl ../genome.final.gff3 geneid_to_mRNAid.txt
cut -f 6 geneid_to_mRNAid.txt >allmRNAID.txt
perl ../scripts/get_mRNA_position.pl allmRNAID.txt ../genome.final.gff3 AT.gff
perl ../scripts/get_fa_by_id.pl allmRNAID.txt ../protein.fa pep.fa
makeblastdb -in pep.fa  -dbtype prot -title pep.fa
blastall -i pep.fa -d pep.fa -e 1e-10  -p blastp  -b 5 -v 5 -m 8 -o AT.blast
MCScanX AT
######################plot##############################################################
mkdir circos
cd circos
perl ../scripts/mcscanX2sync.pl -gff ../MCScanX/AT.gff   -colline ../MCScanX/AT.collinearity    -name genome
samtools faidx  ../genome.fa
awk '{print $1"\t"$2}' ../genome.fa.fai >genome.len
grep -i "^chr" genome.len|sort -k1,1| awk '{print "chr\t-\t"$1"\t"$1"\t0\t"$2"\tchr"NR}'  > karyotype.txt

bedtools  makewindows -w 300000 -g  hap1.genome.len > hap1.genome.window.bed
bedtools  makewindows -w 300000 -g  hap2.genome.len > hap2.genome.window.bed

seqtk subseq ../../01data/hap1.fa  hap1.genome.window.bed  > hap1.genome.window.fa
seqtk subseq ../../01data/hap2.fa  hap2.genome.window.bed  > hap2.genome.window.fa

seqtk comp  hap1.genome.window.fa |awk '{print $1 "\t" ($4+$5)/($3+$4+$5+$6) } ' |awk -F ":|-" '{print $1"\t"$2"\t"$3"\t"$4}'> hap1.gc_density.txt
seqtk comp  hap2.genome.window.fa |awk '{print $1 "\t" ($4+$5)/($3+$4+$5+$6) } ' |awk -F ":|-" '{print $1"\t"$2"\t"$3"\t"$4}'> hap2.gc_density.txt

awk '$3=="gene"{print}'  ../../01data/hap1_genome.gff3 >hap1.gene.gff
awk '$3=="gene"{print}'  ../../01data/hap2_genome.gff3 >hap2.gene.gff

bedtools intersect  -a hap1.genome.window.bed -b hap1.gene.gff   -c -F 0.1  > hap1.genedensity.txt
bedtools intersect  -a hap2.genome.window.bed -b hap2.gene.gff   -c -F 0.1  > hap2.genedensity.txt

awk '$3=="TandemRepeat"{print}'  ../../01data/hap1.all.repeat.gff >hap1.repeat.TandemRepeat.gff
awk '$3=="TandemRepeat"{print}'  ../../01data/hap2.all.repeat.gff >hap2.repeat.TandemRepeat.gff

awk '$3=="Transposon"{print}' ../../01data/hap1.all.repeat.gff >hap1.repeat.TransposonElement.gff
awk '$3=="Transposon"{print}' ../../01data/hap2.all.repeat.gff >hap2.repeat.TransposonElement.gff


grep -i "Class=LTR/Gypsy" ../../01data/hap1.all.repeat.gff >hap1.repeat.LTR_Gypsy.gff
grep -i "Class=LTR/Gypsy" ../../01data/hap2.all.repeat.gff >hap2.repeat.LTR_Gypsy.gff
grep -i "Class=LTR/Copia" ../../01data/hap1.all.repeat.gff >hap1.repeat.LTR_Copia.gff
grep -i "Class=LTR/Copia" ../../01data/hap2.all.repeat.gff >hap2.repeat.LTR_Copia.gff

bedtools coverage -a  hap1.genome.window.bed -b  hap1.repeat.LTR_Gypsy.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > hap1.LTR_Gypsy.coverage.txt
bedtools coverage -a  hap2.genome.window.bed -b  hap2.repeat.LTR_Gypsy.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > hap2.LTR_Gypsy.coverage.txt
bedtools coverage -a  hap1.genome.window.bed -b  hap1.repeat.LTR_Copia.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > hap1.LTR_Copia.coverage.txt
bedtools coverage -a  hap2.genome.window.bed -b  hap2.repeat.LTR_Copia.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > hap2.LTR_Copia.coverage.txt
bedtools coverage -a  hap1.genome.window.bed -b  hap1.repeat.TransposonElement.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > hap1.TransposonElement.coverage.txt
bedtools coverage -a  hap2.genome.window.bed -b  hap2.repeat.TransposonElement.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > hap2.TransposonElement.coverage.txt
bedtools intersect  -a hap1.genome.window.bed -b hap1.repeat.TandemRepeat.gff  -c -F 0.1  > hap1.TandemRepeat.density.txt
bedtools intersect  -a hap2.genome.window.bed -b hap2.repeat.TandemRepeat.gff  -c -F 0.1  > hap2.TandemRepeat.density.txt

cat hap2.TandemRepeat.density.txt hap1.TandemRepeat.density.txt >TandemRepeat.density.txt
cat hap2.TransposonElement.coverage.txt hap1.TransposonElement.coverage.txt >TransposonElement.coverage.txt

cat hap2.LTR_Gypsy.coverage.txt hap1.LTR_Gypsy.coverage.txt > LTR_Gypsy.coverage.txt
cat hap2.genedensity.txt hap1.genedensity.txt > genedensity.txt
cat hap2.gc_density.txt hap2.gc_density.txt > gc_density.txt

####circos--config.txt
echo "
karyotype = hap2_vs_hap1.karyotype.txt #指定染色体文件
chromosomes_units = 500000  #设置长度单位，表示500k长度的序列代表为1u
chromosomes_display_default = yes #默认是将所有的染色体都展示出来
show_tick_labels=yes
show_ticks=yes

##ideogram 染色体和刻度线配置
<ideogram>
    fill=yes
    fill_color=set1-9-qual-1 
    label_font=default
    label_parallel=yes
    label_radius=dims(image,radius)-60p
    label_size=45
    radius=0.90r
    show_label=yes
    <spacing>
        default=10u  #染色体之间间隙
        <pairwise Chr1 chrChr8 >    #设置某两条染色体间隔大一些，方便用于标记label
            spacing = 20u
        </pairwise>
    </spacing>
    stroke_color=dgrey
    stroke_thickness=2p
    thickness=30p
</ideogram>
<ticks> #设置染色体刻度
    color=black
    format=%d
    multiplier=1e-6
    radius=1r
    thickness=2p
    <tick>
        size=10p
        spacing=5u
    </tick>
    <tick>
        color=black
        format=%d
        label_offset=10p
        label_size=25p
        show_label=yes
        size=15p
        spacing=25u
        thickness=4p
    </tick>
</ticks>


############plots##############
track_start = 0.99 #track 起始位置
track_width = 0.06 #track宽度
track_pad   = 0.01 #不同track 之间间隙


<plots>
    <plot>
         type=histogram  # histogram指定绘制直方图 line为折线图，scatter为散点图，heatmap为热图
        fill_color=255,127,80
        color=vlblue
        file=genedensity.txt
        r1  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
        r0  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))

        thickness=2

        <backgrounds> #背景颜色灰色
            <background>
                color = vvlorange
            </background>
        </backgrounds>
    </plot>
   
    <plot>
        type=histogram
        extend_bin = no #直方图是由bin所构成,若bins 在坐标上不相连，最好设置不要将其bins连接到一起。
        file=gc_density.txt
        fill_color=255,140,105
        color=vlblue
        r1  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
        r0  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))
        #stroke_type=both
        thickness=2
        <backgrounds> #背景颜色灰色
            <background>
                color = vlorange
            </background>
        </backgrounds>
    </plot>
    
    <plot>
        type = heatmap
        color   = spectral-5-div
        file = TransposonElement.coverage.txt
        thickness = 5p
        r1  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
        r0  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))
    </plot>
    
    <plot>
        type=histogram
        file=LTR_Copia.coverage.txt
        extend_bin = no #直方图是由bin所构成,若bins 在坐标上不相连，最好设置不要将其bins连接到一起
        fill_color=238,44,44
        color=vlblue
        #r1=0.66r
        #r0=0.56r
        r1  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
        r0  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))

        #stroke_type=both
        thickness=2
        <backgrounds> #背景颜色灰色
            <background>
                color = lorange
            </background>
        </backgrounds>
    </plot>
    
    <plot>
        type=histogram
        file=LTR_Gypsy.coverage.txt
        extend_bin = no #直方图是由bin所构成,若bins 在坐标上不相连，最好设置不要将其bins连接到一起。
        fill_color=238,44,44
        color=vlblue
        r1  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
        r0  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))

        #stroke_type=both
        thickness=2
        <backgrounds> #背景颜色灰色
            <background>
                color = orange
            </background>
        </backgrounds>
    </plot>
   
    <plot>
        type  = heatmap # 指定绘图类型为heatmap
        file  = TandemRepeat.density.txt # 设定数据文件路径
        r1  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))))
        r0  = eval(sprintf(\"%fr\",conf(track_start)-counter(plot)*(conf(track_width)+conf(track_pad))-conf(track_width)))
        #color   = spectral-5-div
        color   = reds-9-seq #reds-9-seq  blues-9-seq
        scale_log_base = 0.5  # 设定 scale_log_base 参数 数据转换 : http://circos.ca/documentation/tutorials/2d_tracks/heat_maps/images
    </plot>


</plots>

<links>  #共线性links
    <link>
        file=hap2_vs_hap1.blocklink.txt
        bezier_radius=0r  # #设置贝塞尔曲线半径，该值设大后曲线扁平。
        bezier_radius_purity=0.75
        color=255,69,0,0.5
        crest=0.6
        #=radius=0.33r. #    手动设置半径
        radius=eval(sprintf(\"%fr\",conf(track_start)-5*(conf(track_width)+conf(track_pad))-conf(track_width)-conf(track_pad))) #自动设置半径，前面6个plot 这里是5，因为索引从0开始
        thickness=2  #设置 link 曲线的粗细
        ribbon=yes  #条带yes,线 no
        #z=20
        <rules>
            <rule> ##染色体内部连线颜色
                color=blue
                condition=var(intrachr)
            </rule>
            <rule>   #染色体间连线颜色
                color=lorange
                condition=var(interchr)
            </rule>

        </rules>
    </link>

</links>

<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
#<<include etc/colors_fonts_patterns.conf>>
#<<include colors.ucsc.conf>>
#<<include colors.hsv.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<image>
<<include etc/image.conf>>
</image>
<<include etc/housekeeping.conf>>

" >config.txt

####circos_polt
circos -conf config.txt -outputdir ./ -outputfile circos
