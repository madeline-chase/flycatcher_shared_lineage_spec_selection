from sys import argv

input = argv[1]

data = []
for line in open(input):
    data.append(line.strip())

filts = ['CG', 'CA', 'TG', 'NG', 'CN']
for trinuc in range(1,len(data),2):
    if data[trinuc][0:2].upper() not in filts and data[trinuc][1:].upper() not in filts:
        print(data[trinuc-1].split(':')[0][1:], int(data[trinuc-1].split(':')[1].split('-')[0])+2,data[trinuc][1],sep = '\t')
