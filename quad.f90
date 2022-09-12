module quad

    !use precisions
    
    implicit none           ! All contained subroutines are also implicit none by default.
    
    private                 ! Default is to make everything private to this module
    
    public :: tetquad  ! tetquad(p,b,wt,npts)  ! order p rules (1 <= p <= 8)
    public :: triquad  ! triquad(p,b,wt,npts)
    public :: edgequad ! edgequad(p,b,wt,npts)
 
    integer, parameter :: RealPrec = 8
!
! Zhang, L., T. Cui, And H. Liu (2009), 
! A Set Of Symmetric Quadrature Rules On Triangles And Tetrahedra,
! Journal Of Computational Mathematics, 27(1), 89–96.
!

    
    contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine tetquad(p,b,wt,npts)
! quadrature rules for tetrahedra
! order p rules (1 <= p <= 8)
! returns: barycentric coordinates of quadrature points,
!          weights, and number of quadratures points
      
      integer :: npts,m,p
      real(RealPrec) :: s,t,w
      real(RealPrec), dimension(*)   :: wt
      real(RealPrec), dimension(4,*) :: b

      select case(p)
      case(1)
         npts=1
         m=0
         s=0.25d0
         w=1.d0
         call S4(m,s,w,b,wt)
      case(2)
         npts=4
         m=0
         s=0.1381966011250105151795413165634361882280d0
         w=0.25d0
         call S31(m,s,w,b,wt)   
      case(3)
         npts=8
         m=0
         s=0.3280546967114266473358058199811974315667d0
         w=0.1385279665118621423236176983756412114150d0
         call S31(m,s,w,b,wt)
         s=0.1069522739329306827717020415706165370512d0
         w=0.1114720334881378576763823016243587885850d0
         call S31(m,s,w,b,wt)
      case(4)
         npts=14
         m=0
         s=0.09273525031089123572902692066289197517846d0
         w=0.07349304311636199595189360213094358276710d0
         call S31(m,s,w,b,wt)
         s=0.3108859192633006945073116831473430813932d0
         w=0.1126879257180159548505940443145759803308d0
         call S31(m,s,w,b,wt)
         s=0.4544962958743503064903641789776029355568d0
         w=0.04254602077708136613167490236965362460141d0
         call S22(m,s,w,b,wt)
      case(5)
         npts=15
         m=0
         s=0.25d0
         w=0.1185185185185185185185185185185185185185d0
         call S4(m,s,w,b,wt)    
         s=0.09197107805272303278884513530051765850491d0
         w=0.07193708377901862001038385490997231731618d0
         call S31(m,s,w,b,wt)     
         s=0.3197936278296299083876254529347764591421d0
         w=0.06906820722627238528062143609531868797482d0
         call S31(m,s,w,b,wt)
         s=0.4436491673103708442589632699891199805416d0
         w=0.05291005291005291005291005291005291005291d0
         call S22(m,s,w,b,wt)
      case(6)
         npts=24
         m=0
         s=0.21460287125915202928883921938628499d0
         w=0.03992275025816749209969062755747998d0
         call S31(m,s,w,b,wt)
         s=0.04067395853461135311557944895641006d0
         w=0.01007721105532064294801323744593686d0
         call S31(m,s,w,b,wt)
         s=0.32233789014227551034399447076249213d0
         w=0.05535718154365472209515327785372602d0
         call S31(m,s,w,b,wt)
         s=0.06366100187501752529923552760572698d0
         t=0.60300566479164914136743113906093969d0
         w=27.d0/560.d0
         call S211(m,s,t,w,b,wt)
      case(7)
         npts=35
         m=0
         s=0.25d0
         w=0.09548528946413084886057843611722638d0
         call S4(m,s,w,b,wt)         
         s=0.31570114977820279942342999959331149d0
         w=0.04232958120996702907628617079854674d0
         call S31(m,s,w,b,wt)
         s=0.05048982259839636876305382298656247d0
         w=0.03189692783285757993427482408294246d0
         call S22(m,s,w,b,wt)
         s=0.18883383102600104773643110385458576d0
         t=0.57517163758700002348324157702230752d0
         w=0.03720713072833462136961556119148112d0
         call S211(m,s,t,w,b,wt)
         s=0.02126547254148324598883610149981994d0
         t=0.81083024109854856111810537984823239d0
         w=0.00811077082990334156610343349109654d0
         call S211(m,s,t,w,b,wt)   
      case(8)
         npts=46
         m=0
         s=0.0396754230703899012650713295393895d0
         w=0.0063971477799023213214514203351730d0
         call S31(m,s,w,b,wt)
         s=0.3144878006980963137841605626971483d0
         w=0.0401904480209661724881611584798178d0
         call S31(m,s,w,b,wt)
         s=0.1019866930627033000000000000000000d0
         w=0.0243079755047703211748691087719226d0
         call S31(m,s,w,b,wt)
         s=0.1842036969491915122759464173489092d0
         w=0.0548588924136974404669241239903914d0
         call S31(m,s,w,b,wt)
         s=0.0634362877545398924051412387018983d0
         w=0.0357196122340991824649509689966176d0
         call S22(m,s,w,b,wt)
         s=0.0216901620677280048026624826249302d0
         t=0.7199319220394659358894349533527348d0
         w=0.0071831906978525394094511052198038d0
         call S211(m,s,t,w,b,wt)
         s=0.2044800806367957142413355748727453d0
         t=0.5805771901288092241753981713906204d0
         w=0.0163721819453191175409381397561191d0
         call S211(m,s,t,w,b,wt)
      case(9:)
         write(*,*) 'Order too large (>8). Using order-1 quadrature.'
         npts=1
         m=0
         s=0.25d0
         w=1.d0
         call S4(m,s,w,b,wt)
      end select

      end subroutine tetquad

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine S4(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(4,*) :: b
      
      n=1
      wt(m+n)=w
      b(:,m+n)=s
      m=m+n
      
      end subroutine S4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine S31(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,t,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(4,*) :: b

      n=4
      t=1.d0 - 3.d0*s
      wt(m+1:m+n)=w
      b(:,m+1)=(/s,s,s,t/)
      b(:,m+2)=(/s,s,t,s/)
      b(:,m+3)=(/s,t,s,s/)
      b(:,m+4)=(/t,s,s,s/)      

      m=m+n
      
      end subroutine S31

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine S22(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,t,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(4,*) :: b

      n=6
      t=0.5d0 - s;
      wt(m+1:m+n)=w      
      b(:,m+1)=(/s,s,t,t/)
      b(:,m+2)=(/s,t,t,s/)
      b(:,m+3)=(/t,t,s,s/)
      b(:,m+4)=(/s,t,s,t/)
      b(:,m+5)=(/t,s,t,s/)
      b(:,m+6)=(/t,s,s,t/)

      m=m+n
      
      end subroutine S22

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine S211(m,s,t,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,t,u,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(4,*) :: b

      n=12
      u=1d0 - 2d0*s-t
      wt(m+1:m+n)=w      
      b(:,m+1)=(/s,s,t,u/)
      b(:,m+2)=(/s,t,s,u/)
      b(:,m+3)=(/t,s,s,u/)
      b(:,m+4)=(/s,s,u,t/)
      b(:,m+5)=(/s,t,u,s/)
      b(:,m+6)=(/t,s,u,s/)
      b(:,m+7)=(/s,u,s,t/)
      b(:,m+8)=(/s,u,t,s/)
      b(:,m+9)=(/t,u,s,s/)
      b(:,m+10)=(/u,s,s,t/)
      b(:,m+11)=(/u,s,t,s/)
      b(:,m+12)=(/u,t,s,s/)
      m=m+n
      
      end subroutine S211

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine S1111(m,s,t,u,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,t,u,v,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(4,*) :: b

      n=24
      v=1d0-s-t-u
      wt(m+1:m+n)=w      
      b(:,m+1)=(/s, t, u, v/)
      b(:,m+2)=(/s, t, v, u/)
      b(:,m+3)=(/s, u, t, v/)
      b(:,m+4)=(/s, u, v, t/)
      b(:,m+5)=(/s, v, t, u/)
      b(:,m+6)=(/s, v, u, t/)
      b(:,m+7)=(/t, s, u, v/)
      b(:,m+8)=(/t, s, v, u/)
      b(:,m+9)=(/t, u, s, v/)
      b(:,m+10)=(/t, u, v, s/)
      b(:,m+11)=(/t, v, s, u/)
      b(:,m+12)=(/t, v, u, s/)
      b(:,m+13)=(/u, s, t, v/)
      b(:,m+14)=(/u, s, v, t/)
      b(:,m+15)=(/u, t, s, v/)
      b(:,m+16)=(/u, t, v, s/)
      b(:,m+17)=(/u, v, s, t/)
      b(:,m+18)=(/u, v, t, s/)
      b(:,m+19)=(/v, s, t, u/)
      b(:,m+20)=(/v, s, u, t/)
      b(:,m+21)=(/v, t, s, u/)
      b(:,m+22)=(/v, t, u, s/)
      b(:,m+23)=(/v, u, s, t/)
      b(:,m+24)=(/v, u, t, s/)
      m=m+n
      
      end subroutine S1111

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine triquad(p,b,wt,npts)
! quadrature rules for triangles
! order p rules (1 <= p <= 10)
! returns: barycentric coordinates of quadrature points,
!          weights, and number of quadratures points
      
      integer :: npts,m,p
      real(RealPrec):: s,t,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(3,*) :: b

      select case(p)
      case(1)
         npts=1
         m=0
         s=1.d0/3.d0
         w=1.
         call R3(m,s,w,b,wt)
      case(2)
         npts=3
         m=0
         s=1.d0/6.d0
         w=1.d0/3.d0
         call R21(m,s,w,b,wt)   
      case(3)
         npts=6
         m=0
         s=0.16288285039589191090016180418490635d0
         w=0.28114980244097964825351432270207695d0
         call R21(m,s,w,b,wt)
         s=0.47791988356756370000000000000000000d0
         w=0.05218353089235368507981901063125638d0
         call R21(m,s,w,b,wt)
      case(4)
         npts=6
         m=0
         s=.44594849091596488631832925388305199d0
         w=.22338158967801146569500700843312280d0
         call R21(m,s,w,b,wt)
         s=.09157621350977074345957146340220151d0
         w=.10995174365532186763832632490021053d0
         call R21(m,s,w,b,wt)
      case(5)
         npts=7
         m=0
         s=1.d0/3.d0
         w=9.d0/40.d0
         call R3(m,s,w,b,wt)    
         s=.10128650732345633880098736191512383d0
         w=.12593918054482715259568394550018133d0
         call R21(m,s,w,b,wt)     
         s=.47014206410511508977044120951344760d0
         w=.13239415278850618073764938783315200d0
         call R21(m,s,w,b,wt)
      case(6)
         npts=12
         m=0
         s=.06308901449150222834033160287081916d0
         w=.05084490637020681692093680910686898d0
         call R21(m,s,w,b,wt)
         s=.24928674517091042129163855310701908d0
         w=.11678627572637936602528961138557944d0
         call R21(m,s,w,b,wt)
         s=.05314504984481694735324967163139815d0
         t=.31035245103378440541660773395655215d0
         w=.08285107561837357519355345642044245d0
         call R111(m,s,t,w,b,wt)
      case(7)
         npts=15
         m=0
         s=.02826392415607634022359600691324002d0
         w=.01353386251566556156682309245259393d0
         call R21(m,s,w,b,wt)
         s=.47431132326722257527522522793181654d0
         w=.07895125443201098137652145029770332d0
         call R21(m,s,w,b,wt)
         s=.24114332584984881025414351267036207d0
         w=.12860792781890607455665553308952344d0
         call R21(m,s,w,b,wt)
         s=.76122274802452380000000000000000000d0
         t=.04627087779880891064092559391702049d0
         w=.05612014428337535791666662874675632d0
         call R111(m,s,t,w,b,wt)   
      case(8)
         npts=16
         m=0
         s=.33333333333333333333333333333333333d0
         w=.14431560767778716825109111048906462d0
         call R3(m,s,w,b,wt)
         s=.17056930775176020662229350149146450d0
         w=.10321737053471825028179155029212903d0
         call R21(m,s,w,b,wt)
         s=.05054722831703097545842355059659895d0
         w=.03245849762319808031092592834178060d0
         call R21(m,s,w,b,wt)
         s=.45929258829272315602881551449416932d0
         w=.09509163426728462479389610438858432d0
         call R21(m,s,w,b,wt)
         s=.26311282963463811342178578628464359d0
         t=.00839477740995760533721383453929445d0
         w=.02723031417443499426484469007390892d0
         call R111(m,s,t,w,b,wt)
      case(9)
         npts=19
         m=0
         s=.33333333333333333333333333333333333d0
         w=.09713579628279883381924198250728863d0
         call R3(m,s,w,b,wt)
         s=.48968251919873762778370692483619280d0
         w=.03133470022713907053685483128720932d0
         call R21(m,s,w,b,wt)
         s=.04472951339445270986510658996627636d0
         w=.02557767565869803126167879855899982d0
         call R21(m,s,w,b,wt)
         s=.43708959149293663726993036443535497d0
         w=.07782754100477427931673935629940396d0
         call R21(m,s,w,b,wt)
         s=.18820353561903273024096128046733557d0
         w=.07964773892721025303289177426404527d0
         call R21(m,s,w,b,wt)
         s=.74119859878449802069007987352342383d0
         t=.22196298916076569567510252769319107d0
         w=.04328353937728937728937728937728938d0
         call R111(m,s,t,w,b,wt)
      case(10)
         npts=25
         m=0
         s=.33333333333333333333333333333333333d0
         w=.08093742879762288025711312381650193d0
         call R3(m,s,w,b,wt)
         s=.42727317884677553809044271751544715d0
         w=.07729858800296312168250698238034344d0
         call R21(m,s,w,b,wt)
         s=.18309922244867502052157438485022004d0
         w=.07845763861237173136809392083439673d0
         call R21(m,s,w,b,wt)
         s=.49043401970113058745397122237684843d0
         w=.01746916799592948691760716329067815d0
         call R21(m,s,w,b,wt)
         s=.01257244555158053273132908502104126d0
         w=.00429237418483282803048040209013191d0
         call R21(m,s,w,b,wt)
         s=.65426866792006614066657009558762790d0
         t=.30804600168524770000000000000000000d0
         w=.03746885821046764297902076548504452d0
         call R111(m,s,t,w,b,wt)
         s=.12280457706855927343012981748128116d0
         t=.03337183373930478624081644177478038d0
         w=.02694935259187995964544947958109671d0
         call R111(m,s,t,w,b,wt)
      case(11:)
         write(*,*) 'Order too large (>10). Using order-1 quadrature.'
         npts=1
         m=0
         s=1.d0/3.d0
         w=1.d0
         call R3(m,s,w,b,wt)
      end select

      end subroutine triquad

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine R3(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(3,*) :: b
      
      n=1
      wt(m+n)=w
      b(:,m+n)=s
      m=m+n
      
      end subroutine R3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine R21(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,t,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(3,*) :: b
      
      n=3
      t=1d0-2d0*s
      wt(m+1:m+n)=w
      b(:,m+1)=(/s,s,t/)
      b(:,m+2)=(/s,t,s/)
      b(:,m+3)=(/t,s,s/)
      m=m+n
      
      end subroutine R21

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     subroutine R111(m,s,t,w,b,wt)

     integer n,m
     real(RealPrec):: s,t,u,w
     real(RealPrec), dimension(*) :: wt
     real(RealPrec), dimension(3,*) :: b

     n=6
     u=1.0d0-s-t
     wt(m+1:m+n)=w
     b(:,m+1)=(/s,t,u/)
     b(:,m+2)=(/t,u,s/)
     b(:,m+3)=(/u,s,t/)
     b(:,m+4)=(/s,u,t/)
     b(:,m+5)=(/t,s,u/)
     b(:,m+6)=(/u,t,s/)
     m=m+n

     end subroutine R111
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine edgequad(p,b,wt,npts)
! gaussian quadrature rules for edges (line segments)
! order p rules (1 <= p <= 10)
! returns: barycentric coordinates of quadrature points,
!          weights, and number of quadratures points
      
      integer :: npts,m,p
      real(RealPrec):: s,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(2,*) :: b

      select case(p)
      case(1)
         npts=1
         m=0
         s=0.5d0
         w=1.d0
         call Q2(m,s,w,b,wt)
      case(2:3)
        npts=2
        m=0
        s=.21132486540518711774542560974902127d0
        w=0.5d0
        call Q11(m,s,w,b,wt)
      case(4:5)
         npts=3
         m=0
         s=0.5d0
         w=4.d0/9.d0
         call Q2(m,s,w,b,wt)
         s=.11270166537925831148207346002176004d0
         w=5.d0/18.d0
         call Q11(m,s,w,b,wt)
      case(6:7)
         npts=4
         m=0
         s=.06943184420297371238802675555359525d0
         w=.17392742256872692868653197461099970d0
         call Q11(m,s,w,b,wt)
         s=.33000947820757186759866712044837766d0
         w=.32607257743127307131346802538900030d0
         call Q11(m,s,w,b,wt)
      case(8:9)
         npts=5
         m=0
         s=0.5d0
         w=128.d0/450.d0
         call Q2(m,s,w,b,wt)    
         s=.04691007703066800360118656085030352d0
         w=.11846344252809454375713202035995868d0
         call Q11(m,s,w,b,wt)     
         s=.23076534494715845448184278964989560d0
         w=.23931433524968323402064575741781910d0
         call Q11(m,s,w,b,wt)
      case(10:11)
         npts=6
         m=0    
         s=.96623475710157601390615077724699730d0
         w=.08566224618958517252014807108636645d0
         call Q11(m,s,w,b,wt)
         s=.83060469323313225683069979750995267d0
         w=.18038078652406930378491675691885806d0
         call Q11(m,s,w,b,wt)     
         s=.61930959304159845431525086084035597d0
         w=.23395696728634552369493517199477550d0
         call Q11(m,s,w,b,wt)
      case(12:13)
         npts=7
         m=0
         s=0.5d0
         w=.20897959183673469387755102040816327d0
         call Q2(m,s,w,b,wt)   
         s=.97455395617137926226309484202392563d0
         w=.06474248308443484663530571633954101d0
         call Q11(m,s,w,b,wt)
         s=.87076559279969721993193238664039420d0
         w=.13985269574463833395073388571188979d0
         call Q11(m,s,w,b,wt)     
         s=.70292257568869858345330320603848073d0
         w=.19091502525255947247518488774448757d0
         call Q11(m,s,w,b,wt)
      case(14:15)
         npts=8
         m=0
         s=.98014492824876811584178043428473650d0
         w=.05061426814518812957626567715498110d0
         call Q11(m,s,w,b,wt)
         s=.89833323870681336979577696823791522d0
         w=.11119051722668723527217799721312044d0
         call Q11(m,s,w,b,wt)
         s=.76276620495816449290886952459462317d0
         w=.15685332293894364366898110099330066d0
         call Q11(m,s,w,b,wt)     
         s=.59171732124782490246973807118009199d0
         w=.18134189168918099148257522463859781d0
         call Q11(m,s,w,b,wt)
      case(16:17)
         npts=9
         m=0
         s=0.5d0
         w=.16511967750062988158226253464348702d0
         call Q2(m,s,w,b,wt)
         s=.98408011975381304491778810145183644d0
         w=.04063719418078720598594607905526183d0
         call Q11(m,s,w,b,wt)
         s=.91801555366331789714971489403486744d0
         w=.09032408034742870202923601562145640d0
         call Q11(m,s,w,b,wt)
         s=.80668571635029519865435101967073709d0
         w=.13030534820146773115937143470931642d0
         call Q11(m,s,w,b,wt)     
         s=.66212671170190446451926900732166830d0
         w=.15617353852000142003431520329222183d0
         call Q11(m,s,w,b,wt)
      case(18:19)
         npts=10
         m=0
         s=.98695326425858586003898200604222603d0
         w=.03333567215434406879678440494666590d0
         call Q11(m,s,w,b,wt)
         s=.93253168334449225536604834421174652d0
         w=.07472567457529029657288816982884867d0
         call Q11(m,s,w,b,wt)
         s=.83970478414951220311716368255743679d0
         w=.10954318125799102199776746711408160d0
         call Q11(m,s,w,b,wt)
         s=.71669769706462359539963297158289208d0
         w=.13463335965499817754561346078473468d0
         call Q11(m,s,w,b,wt)     
         s=.57443716949081560544241300056485999d0
         w=.14776211235737643508694649732566916d0
         call Q11(m,s,w,b,wt)
      case(20:)
         write(*,*) 'Order too large (>20). Using order-1 quadrature.'
         npts=1
         m=0
         s=0.5d0
         w=1.d0
         call Q2(m,s,w,b,wt)
      end select

      end subroutine edgequad

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Q2(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(2,*) :: b
      
      n=1
      wt(m+n)=w
      b(:,m+1)=(/s,s/)
      m=m+n
      
      end subroutine Q2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Q11(m,s,w,b,wt)
      
      integer n,m
      real(RealPrec):: s,t,w
      real(RealPrec), dimension(*) :: wt
      real(RealPrec), dimension(2,*) :: b
      
      n=2
      t=1.d0-s
      wt(m+1:m+n)=w
      b(:,m+1)=(/s,t/)
      b(:,m+2)=(/t,s/)
      m=m+n
      
      end subroutine Q11

end module quad
