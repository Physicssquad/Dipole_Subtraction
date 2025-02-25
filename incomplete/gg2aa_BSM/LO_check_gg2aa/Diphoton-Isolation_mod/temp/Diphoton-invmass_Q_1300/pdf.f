c.----------------------------------------------------------------------
c.                MRST 
c.
c.         x = parton momentum fraction ( 0 < x < 1 )
c.         scale = factorization scale in GeV
c.         iset = 0 for cteq3 lo, 1 for cteq3 nlo msbar
c. output: f filled with parton PROBABILITY distibution functions uses 
c.         Particle Data Group flavor code convention as follows:
c.           tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t
c.           -6  ,-5  ,-4  ,-3  ,-2  ,-1  ,0,1,2,3,4,5,6
c.
c.----------------------------------------------------------------------
      subroutine pdf(x,scale,f)
      implicit double precision (a-h,o-z)
      character name*64
      dimension df(-6:6),f(-6:6)
      common/pdflag/iset
      common/lhapdf/name,mem
c.----------------------------------------------------------------------
      mode=iset
      y=x
      scale1=scale
c      call struct(y,scale1,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
cc      call mrseb(y,scale1,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
c      df(-2)=usea
c      df(-1)=dsea
c      df(+1)=dsea+dnv
c      df(+2)=usea+upv
c      df(-3)=str
c      df(+3)=df(-3)
c      df(-4)=chm
c      df(+4)=df(-4)
c      df(-5)=bot
c      df(+5)=df(-5)
c      df(-6)=0.0d0
c      df(+6)=df(-6)
c      df(0)=glu
cc      endif

       call InitPDFsetByName(name)
       call initpdf(mem)

       call evolvePDF(x,scale,df)

c. Ctq3Pds returns (momentum distribution) = x*(probability distribution)
c. so divide by x to get PROBABILITY distribution
      do 10 i=-6,6
 10   f(i) = df(i)/y

      return
      end
