import os
import re
import matplotlib.pyplot as plt
import optparse

#FileName = raw_input('Please enter your a file name: ')
"""
parse=optparse.OptionParser()
parse.add_option('-g','--gtf',dest='gtf',action='store',metavar='input gtf files',help='enter your transcript gtf)')
parse.add_option('-p','--path',dest='path',action='store',metavar='histone or fantom5 txt file path',help='please enter your txt files')
parse.add_option('-f','--flag',dest='flag',action='store',metavar='',help='input paint tss or histone flag')
(options,args) = parse.parse_args()

#outPutFileName = options.outfile
GTF = options.gtf
PATH = options.path
FLAG = options.flag

"""


def writefile(file_path):
    fhis = open(file_path, 'w')
    ftext = '<=50' + '\t<=100' + '\t<=500' + '\t<=1000' + '\t>1000\n' + \
            '0' + '\t0' + '\t0' + '\t0' + '\t0\n'
    fhis.write(ftext)
    fhis.close()  # 返回gtf中所有的转录本的起始位点

def main(gtf_path,source_path,text_path,pdf_path):                #路径      gtf ，  bed ，  txt
    writefile(text_path)       # 生成统计文件
    transcriptTSS = readGTF(gtf_path)                      #读取GTF文件，晒选得到每个转录本的起始位点的列表
    # print(transcriptTSS)
    stendList = readHistone(transcriptTSS,source_path)     #读histone的文件 返回删选的距离的list
    # stendList = readHG38(transcriptTSS)
    for each in stendList:
        countTSS(each[5],text_path)
    paintTSS(pdf_path,text_path)

def readGTF(genegtfpath):
    transcriptTss = []                      #返回转录本信息（转录起始位点，距离）
    transcriptDistance_list = []              #存储临时转录本信息
    strand = '+'
    chrom = '0'
    num = 0
    with open(genegtfpath, 'r') as f:
        allLines = f.readlines()
        tempName = re.search('transcript_id "(.*?)";', allLines[0]).group(1)
        for eachLine in allLines:
            num +=1                                                                  #记录当前转录的行数
            eachLine_list = eachLine.split('\t')                                              #分割每一行  chr，start，end, strand, geneName,id
            tName = re.search('transcript_id "(.*?)";', eachLine_list[8]).group(1)          #找转录本的名称
            if num == len(allLines):                                                #当为最后一个转录本的时候，将临时名称设空
                tempName = ""
            if tempName != tName:                                                    #设置临时存储的名称
                minData, maxData = countMaxMintranscript(transcriptDistance_list)
                if strand == '+':
                    transcriptTss.append(minData)                           # + 找最小值
                elif strand == '-':
                    transcriptTss.append(maxData)                            # - 找最大值
                transcriptDistance_list.clear()                              #清空上一个的转录本的所有起始距离
                transcriptDistance_list.append((eachLine_list[3], eachLine_list[4], strand, chrom,tName,10001))    #start ,end ,+/-, chr1,transcriptName,distance
                if num == len(allLines):
                    minData, maxData = countMaxMintranscript(transcriptDistance_list)
                    if strand == '+':
                        transcriptTss.append(minData)  # + 找最小值
                    elif strand == '-':
                        transcriptTss.append(maxData)
                tempName = tName
                chrom = eachLine_list[0]
                strand = eachLine_list[6]
            else:
                chrom = eachLine_list[0]
                strand = eachLine_list[6]
                transcriptDistance_list.append((eachLine_list[3], eachLine_list[4], strand, chrom,tName,10001))       #存储每个转录本的信息
    return transcriptTss

def countMaxMintranscript(transcriptDistance_list):               #计算所有的转录本的最大值和最小值
    mintemp = 100000000000
    maxtemp = 0
    for eachtranscript in transcriptDistance_list:                      #获取每个转录本的start & end
        eachxyList =eachtranscript[0]                             #eachxyList[0] =start,   eachxyList[1] = end
        if mintemp - int(eachxyList) >= 0:
            mintemp = int(eachxyList)
            minData = eachtranscript
        if maxtemp < int(eachxyList):
            maxtemp = int(eachxyList)
            maxData = eachtranscript
    return minData,maxData

def readHistone(transcirptList,source_path):
    for each_num in open(source_path, 'r'):                                  #'/source/data/test/Histone_H3K4me1/UCSD.Esophagus.H3K4me1.STL002.bed'
        each_num_list = each_num.split('\t')                                  # ['chr1', '1076372', '1076571', 'HWI-ST216_370:5:1308:8507:9585', '1', '-\n']
        for index,eachTranscirpt in enumerate(transcirptList):                  # 循环遍历存储的转录本的信息        ('130746318', '130746598', '-', '3', 'ENST00000356763')
            # print(index, eachTranscirpt)
            if eachTranscirpt[3] in each_num_list[0]:          # 判断染色体号
                if eachTranscirpt[2] in each_num_list[5]:      # 判断 + / -
                    trandistance = int(eachTranscirpt[0]) + abs(int(eachTranscirpt[1]) -  int(eachTranscirpt[0]))/2       #取外显子的中点与tss中点对比找范围1000以内的值
                    start = each_num_list[1]
                    end = each_num_list[2]
                    if int(start) >= (trandistance - 500) and int(end) <= (trandistance + 500):
                        #找到符合条件的位点，在判断大小
                        standend = abs(int(start) + abs(int(end) - int(start)) / 2 - trandistance)         #histone下的start和end中位置减去转录本的位置，即histone与位点的距离
                        if standend <eachTranscirpt[5]:
                            transcirptList[index] = (eachTranscirpt[0],eachTranscirpt[1],eachTranscirpt[2],eachTranscirpt[3],eachTranscirpt[4],standend)
    return transcirptList

def readHG38(transcirptList,source_path):
    # source_path = os.getcwd() + '/source/' + 'hg38.cage_peak_phase1and2combined_fair_ann.txt.gz.extract.tsv'
    for each_num in open(source_path, 'r'):
        each_num_list = each_num.split('\t')
        chrtext = each_num_list[0]                                                            # chrtext = hg19::chr1:999882..999914,+;hg_64.1
        chr = re.findall('chr(.*?):', chrtext)
        for index, eachTranscirpt in enumerate(transcirptList):        # 循环遍历存储的转录本的信息  ('130746318', '130746598', '-', '3', 'ENST00000356763')
            if eachTranscirpt[3] in chr:  # 判断染色体号
                if eachTranscirpt[2] in chrtext:  # 判断 + / -
                    list1 = re.findall('chr' + eachTranscirpt[3] + ':(.*?),',chrtext)         # list[0] = ['160709055..160709074']
                    each_g38_list = list1[0].split('..')                                        # eachunumlist = ['160709055', '160709074']
                    trandistance = int(eachTranscirpt[0]) + abs(int(eachTranscirpt[1]) -  int(eachTranscirpt[0]))/2       #取外显子的中点与tss中点对比找范围1000以内的值
                    # start = each_g38_list[0]
                    # end = each_g38_list[1]
                    g38distance = int(each_g38_list[0]) + abs(int(each_g38_list[1]) -  int(each_g38_list[0]))/2
                    if g38distance >= (trandistance - 1000) and g38distance <= (trandistance + 1000):
                        # print('11111111')
                        standend = abs(g38distance - trandistance)  # histone下的start和end中位置减去转录本的位置，即histone与位点的距离
                        if standend < eachTranscirpt[5]:
                            transcirptList[index] = (eachTranscirpt[0], eachTranscirpt[1], eachTranscirpt[2], eachTranscirpt[3],eachTranscirpt[4], standend)

    return transcirptList

def countTSS(distanceTss,file_path):         # count distance  -->txt    flag==ture,画小于1000的数
    with open(file_path, 'r') as f1:
        next(f1)
        list1 = f1.readlines()
        line = "".join(list1).split('\t')
        a, b, c, d = line[0], line[1], line[2], line[3]
        e = line[4].replace('\n', '')

    with open(file_path, 'w') as f:
        if distanceTss <= 50:
            num = str(int(a) + 1)
            f.writelines('<=50' + '\t<=100' + '\t<=500' + '\t<=1000' + '\t>1000\n')
            f.writelines(num + '\t' + b + '\t' + c + '\t' + d + '\t' + e + '\n')
        elif distanceTss <= 100:
            num = str(int(b) + 1)
            f.writelines('<=50' + '\t<=100' + '\t<=500' + '\t<=1000' + '\t>1000\n')
            f.writelines(a + '\t' + num + '\t' + c + '\t' + d + '\t' + e + '\n')
        elif distanceTss <= 500:
            num = str(int(c) + 1)
            f.writelines('<=50' + '\t<=100' + '\t<=500' + '\t<=1000' + '\t>1000\n')
            f.writelines(a + '\t' + b + '\t' + num + '\t' + d + '\t' + e + '\n')
        elif distanceTss <= 1000:
            num = str(int(d) + 1)
            f.writelines('<=50' + '\t<=100' + '\t<=500' + '\t<=1000' + '\t>1000\n')
            f.writelines(a + '\t' + b + '\t' + c + '\t' + num + '\t' + e + '\n')
        else:
            num = str(int(e) + 1)
            f.writelines('<=50' + '\t<=100' + '\t<=500' + '\t<=1000' + '\t>1000\n')
            f.writelines(a + '\t' + b + '\t' + c + '\t' + d + '\t' + num + '\n')

def paintTSS(pdf_path,file_path):                                #画tss,histone  的 pdf,txt
    name_list = ['<50', '<100', '<500', '<1000', '>1000']
    with open(file_path, 'r') as f1:
        next(f1)
        list1 = f1.readlines()
        line = "".join(list1).split('\t')
        less50 = int(line[0])
        less100 = int(line[1])
        less500 = int(line[2])
        less1000 = int(line[3])
        more1000 = int(line[4].replace('\n', ''))
        num = less50 + less100 + less500 + less1000 + more1000
    num_list = [less50 / num, less100 / num, less500 / num, less1000 / num, more1000 / num]
    fig1 = plt.figure(1)
    plt.bar(range(len(num_list)), num_list, color='rgb', tick_label=name_list)
    fig1.savefig(pdf_path, dpi=150)
    plt.close()

def chooseFile(gtf_path,source_path,flag):
    if flag == 'fantom5':          #if flag == True,用户产生fantom5 的数据
        his_files = os.listdir(source_path)  # 遍历histone列表获取文件名
        for fi in his_files:  # 遍历得到每一个文件名
            tss_path = 'My result/' + fi.split('.bed')[0]
            writefile(tss_path + '.txt')
            # 返回转录本信息列表
            tsstranscriptList = readGTF(gtf_path)
            tssList = readHistone(tsstranscriptList, source_path+ '/' + fi)
            for each in tssList:
                countTSS(each[5], tss_path + '.txt')
            paintTSS(tss_path + '.pdf', tss_path + '.txt')
    elif flag == 'histom':
        his_files = os.listdir(source_path)  # 遍历histone列表获取文件名
        for fi in his_files:  # 遍历得到每一个文件名
            # print(fi.split('.bed')[0])
            H3K4_path = 'My result/' + fi.split('.bed')[0]   # fi.split('.bed')[0] = 文件名
            writefile(H3K4_path + '.txt')
            H3K4transcriptList = readGTF(gtf_path)  # 返回转录本信息列表
            histoneH3K4 = readHistone(H3K4transcriptList, source_path + '/' + fi)  # 更新转录本信息，添加转录本距离
            for each in histoneH3K4:
                countTSS(each[5], H3K4_path + '.txt')
            paintTSS(H3K4_path + '.pdf',H3K4_path + '.txt')


if __name__ == '__main__':
    chooseFile(GTF,PATH,FLAG)
    # chooseFile('source/K510_test.gtf','histom5','fantom5')
