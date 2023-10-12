import sympy as smp




def Hk2html(ParaIn,HkMatSp,folderName="./"):
    """
    :param ParaIn
    :param HkMatSp: sympy matrix of Hk
    :param folderName: output folder of HTML file
    :return: print Hk to html file
    """

    #write entire matrix
    nRow,nCol=HkMatSp.shape
    HkReIm=smp.Matrix.zeros(nRow,nCol)
    for i in range(0,nRow):
        for j in range(0,nCol):
            reTmp,imTmp=HkMatSp[i,j].rewrite(smp.cos).expand().as_real_imag()
            HkReIm[i,j]=reTmp+smp.I*imTmp

    content="""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">"""\
            +"\n"\
            +"<html>" \
            + "\n" \
            +"<head>" \
            + "\n" \
            +"<title>Mathedemo</title>" \
            + "\n" \
            +"""<script type="text/x-mathjax-config">""" \
            + "\n" \
            +"""  MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});""" \
            + "\n" \
            +"</script>" \
            + "\n" \
            +"""<script type="text/javascript" """ \
            + "\n" \
            +"""  src="http://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS-MML_HTMLorMML">""" \
            + "\n" \
            +"</script>" \
            + "\n" \
            +"</head>" \
            + "\n" \
            +" " \
            + "\n" \
            +"<body>" \
            + "\n" \
            +"<h2>"+"Hamiltonian of "+ParaIn["Name"]+" in momentum space"+"</h2>" \
            + "\n" \
            +"$"+smp.latex(HkReIm)+"$"

    fPtr = open(folderName+"/" + ParaIn["Name"] + "_Hk.html","w+")
    fPtr.write(content)
    fPtr.write("   \n")
    fPtr.write("<h2>Display each element</h2>")
    for i in range(0,nRow):
        for j in range(0,nCol):
            fPtr.write("Matrix element ["+str(i)+", "+str(j)+"] is $"+smp.latex(HkReIm[i,j])+"$\n")
            fPtr.write("  \n")
    fPtr.close()