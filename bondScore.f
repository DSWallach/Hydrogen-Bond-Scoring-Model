      program bondScore

c.....This program assigns a hydrogen bond score to a DNA sequence
c     to produce the score this program requires the dist.dat file
c     produced by running cpptraj_dist on the desired DNA sequence
c     04.08.15 David Selwyn Wallach, Kelly M. Thayer

      
      integer until             ! Used to create a repear-until loop
      integer max               ! Also used for looping
      integer until2            !
      integer max2              !
      integer ID(9)             ! Bond input identites, Integer form
      real length(9)            ! readin sequence
      character*8 resid         ! string "Residue "
      character*8 in            ! for input
      character*8 lengths       ! string " Lengths"
      character*8 scores        ! string "  Scores"
      character*4 filename      ! holds the file name w/ no data type
      character*9 bond(9)       ! Bond names
      character*9 knownID(9)    ! Known bod input identities
      character*1 iden(9)       ! Bond input identities
      character*13 outfile      ! output file name
      character*12 infile       ! input file name
      character*43 sumString    ! used to output the sm of all the scores
      
c.....Bins
      real bin1                 ! < 2.635
      real bin2                 ! < 2.7
      real bin3                 ! < 2.765
      real bin4                 ! < 2.81
      real bin5                 ! < 3.09
      real bin6                 ! < 3.2
      real bin7                 ! < 3.3
      real bin8                 ! < 3.4
      real bin9                 ! < 3.5

c.....Score
      real score(9)             ! The score given to the input file
      real bin1score            ! A interger representing the percentage 
      real bin2score            ! of 1 that will be added to score
      real bin3score            ! for each bin the bond could be contained 
      real bin4score            ! in
      real bin5score            !
      real bin6score            !
      real bin7score            !
      real bin8score            !
      real bin9score            !
      real sumScore             ! holds the sum of the scores of each bond

c.....Initialization of Bins
 98   bin1=2.635
      bin2=2.7
      bin3=2.765
      bin4=2.81
      bin5=3.09
      bin6=3.2
      bin7=3.3
      bin8=3.4
      bin9=3.5

c.....Intialization of  bin#score
 97   bin1score=1.98019801980198
      bin2score=7.60869565217391
      bin3score=15.6626506024096
      bin4score=21.6867469879518
      bin5score=25
      bin6score=20.6896551724138
      bin7score=13.7931034482759
      bin8score=7.77777777777778
      bin9score=3.33333333333333

c.....Initialization of bond(9)
 96   bond(1)= "A Lys 120"
      bond(2)= "B Lys 120"
      bond(3)= "  Ser 241"
      bond(4)= "  Ala 276"
      bond(5)= "  Cys 277"
      bond(6)= "A Arg 280"
      bond(7)= "B Arg 280"
      bond(8)= "  Arg 283"
      bond(9)= "  Arg 273"
      resid= "Residue:"

c.....Intialization of length(9), necessary for when the input file
c      has fewer than nine bonds 
 95   length(1)=0
      length(2)=0
      length(3)=0
      length(4)=0
      length(5)=0
      length(6)=0
      length(7)=0
      length(8)=0
      length(9)=0

c.....Initializing iden(9) and others
 94   iden(1)="9"
      iden(2)="9"
      iden(3)="9"
      iden(4)="9"
      iden(5)="9"
      iden(6)="9"
      iden(7)="9"
      iden(8)="9"
      iden(9)="9"
      until=1
      max=9
      until2=1
      max2=9

      print*, "Welcome to bondScore!"      
      print*, "Please input file name (xxx_dist.dat)"
      read*, infile
      print*, "Your input file is ", infile
      
c.....Source input
      open(unit=19,file=infile,status='old')
      read(19,503)label,iden(1),iden(2),iden(3),iden(4),iden(5),
     1     iden(6),iden(7),iden(8),iden(9)
      read(19,501)in,length(1),length(2),length(3),
     1     length(4),length(5),length(6),length(7),length(8),length(9)
      
c.....Converting character inden(9) to integer ID(9)
      print*,"Working before convert"
 150  CONTINUE
      read(iden(until),*)ID(until)
      until=until+1
      if((ID(until)==until-1).AND.until.LT.9)then
         goto 150
      elseif(until.LT.9)then
         until2=until
         goto 152
      else 
         goto 160
      endif

c.....Fixing missing input
 152  DO WHILE(.NOT.(10 - until2==1+ID(10 - until2))) 
         iden(10 - until2)=iden((10 - until2) - 1)
         length(10 - until2)=length((10 - until2) - 1)
         ID(10 - until2)=ID((10 - until2) - 1)
         until2=until2+1
      END DO
      goto 150
      
 160  CONTINUE
c.....do loop for binning each length
      print*,"Working before DO loop"
      DO 50, i = 1, 9
         print*, "working before first if",ID(i),i
         if (.NOT.(ID(i+1)==9.AND.8==i))then
            print*, "Working after first if"
            if (length(i).LT.bin1)then
               score(i)=bin1score
            elseif(length(i).LT.bin2)then
               score(i)=bin2score
            elseif(length(i).LT.bin3)then
               score(i)=bin3score
            elseif(length(i).LT.bin4)then
               score(i)=bin4score
            elseif(length(i).LT.bin5)then
               score(i)=bin5score
            elseif(length(i).LT.bin6)then
               score(i)=bin6score
            elseif(length(i).LT.bin7)then 
               score(i)=bin7score
            elseif(length(i).LT.bin8)then
               score(i)=bin8score
            elseif(length(i).LT.bin9)then
               score(i)=bin9score
            else
               score(i)=0
               print*, "At", bond(i),"Length is outside of bin range"
            endif
         else
            print*, "Working after first else"
            score(i)=0
            print*, "At", bond(i),"No hydrogen bond was formed"
c     goto 93
         endif
 50   CONTINUE

c......Incase a bond was not formed
c     print*,"Working before convert + 1"
c     93   read(iden(i+1),*)ID(i)
      
c.....Sumation of scores
      print*, "Working after DO loop"
 89   sumScore=score(1)+(score(2)+(score(3)+(score(4)+(score(5)+
     1     (score(6)+(score(7)+(score(8)+score(9))))))))
      sumString="The sum of the scores of all the bonds is"
      filename=infile
      outfile=filename//"score.out"
      lengths=" Lengths"
      scores= " Scores "

      print*, "Your output file is ", outfile      
c.....Data output
      open(unit=20,file=outfile)
      write(20,500)resid,bond(1),bond(2),bond(3),bond(4),bond(5),
     1     bond(6),bond(7),bond(8),bond(9)
      write(20,501)lengths,length(1),length(2),length(3),
     1     length(4),length(5),length(6),length(7),length(8),length(9)
      write(20,501)scores,score(1),score(2),score(3),score(4),
     1     score(5),score(6),score(7),score(8),score(9)
      write(20,502),sumString,sumScore

 500  format(A8,9(1x,A12))
 501  format(A8,9(6x,f7.4))
 502  format(A43,f9.4)
 503  format(A8,9(12x,A1))
      end
