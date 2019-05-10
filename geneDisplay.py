import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import re
import optparse
import paint_H3K4

"""
parse=optparse.OptionParser()
parse.add_option('-g','--gtf',dest='gtf',action='store',metavar='input gtf files',help='enter your transcript gtf)')
parse.add_option('-q','--quant',dest='txt',action='store',metavar='txt file name',help='please enter your txt files')
parse.add_option('-i','--id',dest='id',action='store',metavar='',help='input gene ensg id')
parse.add_option('-p','--path',dest='path',action='store',metavar='input histone listdir',help='please enter your histone listdir')

(options,args) = parse.parse_args()
"""

class showGene(object):
    def __init__(self, gtfText_path, Gene_name, idtext_path,histonePath):
        self.genename = ''
        self.gene_gtf = ''
        self.transcript_num = 0
		# 某基因下的转录本数目
        self.ucsctransNum = 0
		#UCSC 转录本数目
        self.gene_gtf = gtfText_path
		#input GTF
        self.Title = self.transferName(Gene_name.upper())
        self.histone_path = histonePath
		#输入的histone文件夹的名称
        self.tsssource_path = os.path.split(os.path.realpath(__file__))[0]+'/'+'source/hg38.cage_peak_phase1and2combined_fair_ann.txt.gz.extract.tsv'
        try:
            self.genename = self.Title[2]
			#user input Gene_name
        except:
            print("genename no exist,please check and retry")
            self.genename = ''
        self.show(idtext_path)

    def Read_WriteFile(self,dataGTF_path,genetxt_path,flag):
        transcript_list = []
        fp = open(genetxt_path, 'w')
        for eachLine in open(dataGTF_path, 'r'):
            eachLine_list = eachLine.split('\t')  # split everyline
            name = re.search('gene_id "(.*?)";', eachLine_list[-1]).group(1)  # search everyline's Gene_name in GTF files
            if self.genename in name:  # search the name is equal of  input
                fp.write(eachLine)  # find the date input on the rsult of geneTXT file
                if eachLine_list[2] == 'exon':
                    transcript_name = re.search('transcript_id "(.*?)";', eachLine_list[-1]).group(1)  # search transcript_name under the genename
                    transcript_name = transcript_name.split('.')
                    if transcript_name[0] not in transcript_list:
                        transcript_list.append(transcript_name[0])
                        if flag == False:
                            self.ucsctransNum +=1
                        else:
                            self.transcript_num += 1
        fp.close()
        return transcript_list

    def show(self, idtext_path):
        global genetxt_path,ucsc_path
        if self.gene_gtf == '' or self.genename == '':
            return
        result_path = self.genename         #genenameTXTfile include pdf and gene file
        if not os.path.exists(result_path):
		#if don't have result path ,create file
            os.makedirs(result_path)
        genetxt_path = result_path + '/' + self.genename + '.txt';#gene file path
        ucsc_path = result_path + '/UCSC' + self.genename + '.txt'
        ucsc_filepath = os.path.split(os.path.realpath(__file__))[0]+'/'+'source/' + 'gtfAnnotation.gtf'
        transcript_list = self.Read_WriteFile(self.gene_gtf,genetxt_path,True);#store this specile gene's transcript_list
        # print(transcript_list)
        ucsctranscrpt_list = self.Read_WriteFile(ucsc_filepath,ucsc_path,False);#store ucsc list
        if transcript_list == []:
            print("genename no exist,please check and retry")
        else:
            lista = [];#recorde ax1 label
            for eachLine in ucsctranscrpt_list:
                lista.append(' ')
            for eachnum in transcript_list:
                if idtext_path != '':
                    trLabel = '(' + str(self.read_transcriptNum(idtext_path, eachnum)) + ')'
                else:
                    trLabel = ''
                lista.append(trLabel)
            self.paint(result_path,transcript_list,ucsctranscrpt_list,lista)

    def paint(self,result_path,transcript_list,ucsctranscrpt_list,lista):
        line_width = 5;# the line width of picture
        frontsize = 6
        quver_width = 0.003
        if self.transcript_num >= 21:
		#if transcript_num>21,line_width = 4
            line_width = 4
            frontsize = 4
            quver_width = 0.005
        elif self.transcript_num + self.ucsctransNum >=40:
            line_width = 3
            frontsize = 3
            quver_width = 0.006
        fig = plt.figure(1)
        num = 0;#所有转录本的数目
        traNum = -1;#所有转录本start&end的数目
        pdf_path = result_path + '/' + self.genename + '.pdf'
        temporaryFile = result_path + '/temp.txt'; #temparay file to store src & f  used it destory
        src = open(genetxt_path, "r")
        f = open(ucsc_path, 'r')
        temp = open(temporaryFile, "w")
        temp.write(f.read())
        temp.write(src.read())
        f.close()
        src.close()
        temp.close()
        trStart = []
        trEnd = []; #count transcript min start max end
        for eachLine in  open(temporaryFile, 'r'):                  #first : searching the total data of gene transcript for max&min transcript
            eachLine_list = eachLine.split('\t')
            trStart.append(eachLine_list[3])
            trEnd.append(eachLine_list[4])
        transcript_start = min(trStart)
        transcript_end = max(trEnd)
        trStart.clear()
        trEnd.clear()
        fp = open(temporaryFile, 'r')
        allLines = fp.readlines()
        for eachLine in allLines:
            traNum = traNum + 1
            eachLine_list = eachLine.split('\t')
            trStart.append(eachLine_list[3])
            trEnd.append(eachLine_list[4])
            tName = re.search('transcript_id "(.*?)";', allLines[traNum]).group(1)
            if (traNum + 1) < len(allLines):
                tName1 = re.search('transcript_id "(.*?)";', allLines[traNum + 1]).group(1)
            else:
                tName1 = ""
            if eachLine_list[2] == 'exon':
                if eachLine_list[6] == '+':
                    arr = '4'
                else:
                    arr = '3'
                #左纵坐标轴ax
                ax = fig.add_axes([0.2, 0.2, 0.5, 0.6])
                ax.set_xlim(int(transcript_start) - 10, int(transcript_end) + 10)
                ax.set_ylim(-0.5, self.transcript_num + self.ucsctransNum +2)
                ax.set_xticks(np.linspace(int(transcript_start) - 1000, int(transcript_end) + 1000, 2))
                ax.set_yticks([0.1] + list(range(1, self.transcript_num + self.ucsctransNum + 2)))
                ax.set_yticklabels(['HISTONE'] + ['FANTOM5'] + ucsctranscrpt_list + transcript_list)   #y轴标签
                ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace
                (int(transcript_start) - 10, int(transcript_end) + 10, 2)]])
                ax.spines['top'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.get_xaxis().set_tick_params(direction='out')
                ax.tick_params(axis=u'y', which=u'both', length=0)
                #copy 左纵坐标，得到相同的右坐标
                ax1 = ax.twinx()
                ax1.set_ylim(-0.5, self.transcript_num+ self.ucsctransNum + 2)
                ax1.set_yticks([0.1] + list(range(1, self.transcript_num+ self.ucsctransNum + 2)))
                ax1.set_yticklabels([' '] + [' '] + lista)
                ax1.spines['top'].set_visible(False)
                ax1.spines['left'].set_visible(False)
                ax1.spines['right'].set_visible(False)
                ax1.tick_params(axis=u'y', which=u'both', length=0)
                #设置字体大小
                for label in ax1.get_yticklabels():
                    label.set_fontsize(frontsize)
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(frontsize)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(frontsize)
                if num == 0:
                    line1 = [int(transcript_start) - 1000, num + 0.1], [int(transcript_end) + 1000, num + 0.1]
                    (line1_xs, line1_ys) = zip(*line1)
                    ax.add_line(lines.Line2D(line1_xs, line1_ys, linewidth=0.5, c='black'))
                    gtf_path = self.genename + '/' + self.genename + '.txt';#gtf截取的基因文件名
                    his_files = os.listdir(self.histone_path);#遍历histone列表获取文件名
                    colorList = ['red','blue','green','orange','pink','black'];#存储颜色
                    colornum = 0
                    for fi in his_files:                             #遍历得到每一个文件名，绘制箭头
                        H3K4_path = self.genename + '/'+fi.split('.bed')[0]+'.txt';
                        paint_H3K4.writefile(H3K4_path);
                        H3K4transcriptList = paint_H3K4.readGTF(gtf_path);# 返回转录本信息列表
                        histoneH3K4 = paint_H3K4.readHistone(H3K4transcriptList, self.histone_path + '/' + fi);#更新转录本信息，添加转录本距离
                        for each in histoneH3K4:
                            paint_H3K4.countTSS(each[5], H3K4_path)
                        for eachHis in histoneH3K4:  # 便遍历列表，绘制箭头
                            if eachHis[5] < 1000:
                                int(eachHis[0]) + eachHis[
                                    5];# ('70301925', '70302067', '-', '2', 'ENST00000460307', 46.5)
                                line = [(int(eachHis[0]) + eachHis[5] - 1, num + 0.1),
                                        (int(eachHis[1]) - eachHis[5] + 1, num + 0.1)]
                                (line_xs, line_ys) = zip(*line)
                                ax.quiver(line_xs, line_ys, [0, 0], [0, 1], color= colorList[colornum] , width= quver_width,scale=25)  # 绘制箭头
                        colornum +=1
                    num = num + 1      #行数+1
                if num == 1:                     #绘制tss
                    line1 = [int(transcript_start)-1000, num + 0.1], [int(transcript_end)+1000, num + 0.1]
                    (line1_xs, line1_ys) = zip(*line1)
                    ax.add_line(lines.Line2D(line1_xs, line1_ys, linewidth=0.5, c='black'))
                    tss_path = self.genename + '/tss.txt'
                    gtf_path = self.genename + '/' + self.genename + '.txt'
                    paint_H3K4.writefile(tss_path);
                    # 返回转录本信息列表
                    tsstranscriptList = paint_H3K4.readGTF(gtf_path)
                    tssList = paint_H3K4.readHistone(tsstranscriptList, self.tsssource_path)
                    for each in tssList:
                        paint_H3K4.countTSS(each[5], tss_path)
                    for eachTss in tssList:
                        if eachTss[5] < 1000:
                            int(eachTss[0]) + eachTss[5]  # ('70301925', '70302067', '-', '2', 'ENST00000460307', 46.5)
                            line = [(int(eachTss[0]) + eachTss[5] - 1, num + 0.1),
                                    (int(eachTss[1]) - eachTss[5] + 1, num + 0.1)]
                            (line_xs, line_ys) = zip(*line)
                            ax.quiver(line_xs, line_ys, [0, 0], [0, 1], color='red', width=quver_width, scale=25)
                    num = num + 1
                if num-1 <= self.ucsctransNum :
                    setcolor = '#0013B3'
                else:
                    setcolor = '#9F0131'
                line1 = [(int(eachLine_list[3]), num), (int(eachLine_list[4]), num)]
                (line1_xs, line1_ys) = zip(*line1)
                ax.add_line(lines.Line2D(line1_xs, line1_ys,
                                         solid_capstyle='butt', solid_joinstyle='miter',
                                         linewidth=int(line_width), alpha=1,
                                         antialiased=False, color=setcolor))
                if tName != tName1:
                    minLine = int(min(trStart))
                    maxLine = int(max(trEnd))
                    trline = [(minLine, num),
                              (maxLine, num)]
                    (trline_xs, trline_ys) = zip(*trline)
                    ax.add_line(lines.Line2D(trline_xs, trline_ys,
                                             solid_capstyle='butt', solid_joinstyle='miter',
                                             linewidth=0.5, alpha=1,
                                             antialiased=False, color=setcolor))
                    #画横轴上的方向，设定为30个箭头->-->--->
                    addtext = (int(transcript_end) - int(transcript_start)) / 30
                    while (minLine < maxLine):
                        ax.scatter(minLine, num, marker=arr, linewidth=0.1, c='blue', s=15, alpha=0.5)
                        minLine += addtext
                    num = num + 1
                    trStart.clear()
                    trEnd.clear()

            fig.suptitle('\n\n\nchr' + str(eachLine_list[0]) + ': ' + self.genename + '\n Entrez Gene:' +
                         self.Title[0] +'/Symbol ID:'+ self.Title[1], fontsize=10)
        fig.savefig(pdf_path, dpi=150)
        plt.close()
        self.paintTSS(self.genename + '/tss.pdf',tss_path)
        his_files = os.listdir(self.histone_path);# 遍历histone列表获取文件名  绘制柱状图
        for fi in his_files:
            self.paintTSS(self.genename +  '/' + fi.split('.bed')[0] + '.pdf', self.genename + '/' + fi.split('.bed')[0] + '.txt')

    def read_transcriptNum(self,text_path,geneName):             #读取id2id的文件，绘制右纵坐标标签
        total = 0;
        trNumList = []
        for each_num in open(text_path, 'r'):
            if not each_num.startswith('#'):
                each_num_list = each_num.split('\t')
                trname = each_num_list[1]
                if trname == geneName:
                    trNumList.append(each_num_list[0])
        for each_num in trNumList:
            numlist = re.findall(r"\d+",re.search('/f(.*?)p',each_num).group())
            num = int(numlist[0])
            total = total + num
        return total

    def transferName(self,Gene_name):        #用户输入基因名，遍历name文件找该基因的别名
        Title = []
        genename_path = os.path.split(os.path.realpath(__file__))[0]+'/'+'source/' + 'name.txt'
        f = open(genename_path, 'r')
        allLines = f.readlines()
        for eachLine in allLines:
            eachLine_list = eachLine.split('\t')
            if Gene_name == eachLine_list[0]:
                Title = eachLine_list
            elif Gene_name == eachLine_list[1]:
                Title = eachLine_list
            elif Gene_name == eachLine_list[2].replace('\n', ''):
                Title = eachLine_list
        Title[2] = Title[2].replace('\n', '')
        print(Title)
        f.close()
        return Title

    def paintTSS(self,pdf_path,txtpath):                                #画tss,histone  的 pdf,txt
        name_list = ['<50', '<100', '<500', '<1000', '>1000']
        with open(txtpath, 'r') as f1:
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

"""
def main():
    #基因文件gtf，   基因名 ， 右坐标轴标签 ， HH3K4m1_path ， H3K4m3_path  ，H3K27ac_path
    test = showGene('source/K510_3rd.gtf', 'ezr', 'source/K510_3rd.id2id.txt','histone')

if __name__ == '__main__':
    main()

"""
def geneDisplay(GTF, GE_ID, TXT, HIS_PATH):
	showGene(GTF,GE_ID,TXT,HIS_PATH);




