def copyLines(oldName, maxLines, newName):
    arc = open(oldName)
    lines = arc.readlines()
    linesNew = []
    for i in range(0, min(len(lines), maxLines)):
        linesNew = linesNew + [lines[i]]
    arc.close()
    arc = open(newName, 'w')
    arc.writelines(linesNew)
    arc.close()

copyLines('ShuTransfer.py', 500, 'ShuCopy.py')
