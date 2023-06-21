!! 二次元Helmholtz方程式の外部Neumann問題を境界要素法で解くプログラム
program main
  implicit none

  ! 数学定数
  complex(8),parameter :: ione=(0.d0,1.d0) !< 虚数単位

  ! Gauss-Legendre 数値積分
  integer,parameter :: ng=10 !< ガウス積分の分点の数
  real(8) :: gzai(ng) !< 分点
  real(8) :: wei(ng) !< 重み

  ! Helmholtz 方程式
  real(8) :: wn=3.d0 !< 波数(wavenumber)=k
  real(8),parameter :: xs(2)=[-1.5d0, 0.d0] !点源の位置
  
  ! 境界要素メッシュ
  integer :: n !< 境界要素数
  real(8),allocatable :: y(:,:) !< y(i,j)=j番節点のx_i座標
  real(8),allocatable :: x(:,:) !< x(i,j)=j番要素の選点(中点)のx_i座標
  integer,allocatable :: nd(:,:) !< nd(1,j): j番要素の始点の節点番号, nd(2,j): j番要素の終点の節点番号
  real(8),allocatable :: an(:,:) !< an(i,j): j番要素の法線ベクトル
  real(8),allocatable :: le(:) !< le(j): j番要素の長さ

  ! 境界要素法
  complex(8),allocatable :: hmat(:,:) !< 境界要素法に現れる係数行列
  complex(8),allocatable :: p(:) !< i番選点における複素音圧(境界値)
  real(8) :: yxi(2) !< 式(86)のy

  ! 内点計算
  real(8) :: amin(2), amax(2) !< 内点を [amin(1),amax(1)]x[amin(2),amax(2)] に配置
  integer :: nip(0:2) ! nip(1)=x1方向の内点の数-1, nip(2)=x2方向の内点の数-1, nip(0)=内点の総数
  real(8),allocatable :: xip(:,:) !< xip(i,j): j番内点のx_i座標
  complex(8),allocatable :: pip(:) !< pip(i): i番内点の複素音圧
  
  ! その他の変数
  integer :: i, j, k, itmp, ig, ip
  real(8) :: xt, yt, r, dotp
  complex(8) :: hank, tmp
  real(8) :: intensity
  
  !---------------------------------------------------------------------------
  ! Gauss-Legendre数値積分の分点と重みを設定
  gzai=(/-0.97390652851717174d0,&
       -0.86506336668898454d0,&
       -0.67940956829902444d0,&
       -0.43339539412924721d0,&
       -0.14887433898163122d0,&
       0.14887433898163122d0,&
       0.43339539412924721d0,&
       0.67940956829902444d0,&
       0.86506336668898454d0,&
       0.97390652851717174d0/)

  wei=(/6.6671344308688138d-2,&
       0.14945134915058059d0,&
       0.21908636251598204d0,&
       0.26926671930999635d0,&
       0.29552422471475287d0,&
       0.29552422471475287d0,&
       0.26926671930999635d0,&
       0.21908636251598204d0,&
       0.14945134915058059d0,&
       6.6671344308688138d-2/)

  !---------------------------------------------------------------------------
  ! 境界要素メッシュを読み込み
  open(1,file="mesh.txt")
  read(1,*) n !節点数
  allocate(y(2,n)) !節点数が決まったので節点座標を格納する配列を確保
  do i=1,n
     read(1,*) itmp, y(1,itmp), y(2,itmp) !順番に節点座標を読み込み
  end do
  read(1,*) n !要素数
  allocate(nd(2,n)) !要素数が決まったので要素の始点と終点の節点番号を格納する配列を確保
  do i=1,n
     read(1,*) itmp, nd(1,itmp), nd(2,itmp) !順番に節点座標を読み込み
  end do
  close(1)

  !---------------------------------------------------------------------------
  ! 選点（要素中心）座標を計算
  allocate(x(2,n))
  do i=1,n
     x(:,i)=(y(:,nd(1,i))+y(:,nd(2,i)))*0.5d0
  end do

  !---------------------------------------------------------------------------
  ! 法線ベクトルと要素長さを計算
  allocate(an(2,n)) !法線ベクトル
  allocate(le(n)) !要素長さ
  do i=1,n
     xt=y(1,nd(2,i))-y(1,nd(1,i)) !要素に沿うベクトル(xt, yt)を
     yt=y(2,nd(2,i))-y(2,nd(1,i)) ! 計算
     le(i)=????? !i番要素の長さ
     xt=xt/le(i) !要素に沿う単位ベクトルのx1成分
     yt=yt/le(i) !要素に沿う単位ベクトルのx2成分
     an(1,i)=yt !法線ベクトルは要素に沿うベクトル
     an(2,i)=-xt !に直交するように設定すれば良い（向きに注意）
  end do
  
  !---------------------------------------------------------------------------  
  ! 内点の数を設定
  nip(1)=99 ! x1方向の内点の数-1
  nip(2)=99 ! x2方向の内点の数-1
  nip(0)=(nip(1)+1)*(nip(2)+1) !内点の総数

  ! 内点の存在する範囲を設定
  amin(1)=-2.d0
  amax(1)=+2.d0
  amin(2)=-2.d0
  amax(2)=+2.d0
  
  ! 内点
  allocate(xip(2,nip(0)))
  nip(0)=0
  do j=0,nip(2)
     yt=(amax(2)-amin(2))*dble(j)/dble(nip(2))+amin(2)
     do i=0,nip(1)
        xt=(amax(1)-amin(1))*dble(i)/dble(nip(1))+amin(1)
        nip(0)=nip(0)+1
        xip(1,nip(0))=xt
        xip(2,nip(0))=yt
     end do
  end do

  !---------------------------------------------------------------------------
  ! 境界要素法の行列を計算(86)式
  allocate(hmat(n,n))
  hmat(:,:)=cmplx(0.d0,0.d0) !初期化
  do i=1,n !選点のループ
     do j=1,n !要素のループ
        do ig=1,ng !ガウス積分のループ
           yxi(:)=????? !式(86)
           r=????? ! r=|x-y(ξ)| 式(87)
           dotp=????? ! dotp=(x-y(ξ))・n 式(87)
           hank=????? ! hank=H_1^(1)(k|x-y(ξ)|) 式(87)
           hmat(i,j)=hmat(i,j)+0.125d0*ione*wn*le(j)*hank*dotp/r*wei(ig) !(87)式
        end do
     end do
  end do
  do i=1,n
     hmat(i,i)=hmat(i,i)+0.5d0 !係数行列に式(85)の1/2δijを足す
  end do

  !---------------------------------------------------------------------------
  ! 境界選点上の入射波を作成しpに格納
  allocate(p(n))
  do i=1,n
     r=sqrt(dot_product(x(:,i)-xs(:),x(:,i)-xs(:))) ! 点源と境界選点の距離
     p(i)=0.25d0*ione*(bessel_j0(wn*r)+ione*bessel_y0(wn*r)) ! 式(85)の右辺をp(i)に代入
  end do

  !---------------------------------------------------------------------------
  ! 連立一次方程式(85)を解いて境界選点上のpを求める
  do k=1,n-1 !前進消去
     do i=k+1,n 
        tmp=hmat(i,k)/hmat(k,k)
        do j=k+1,n
           hmat(i,j)=hmat(i,j)-tmp*hmat(k,j)
        end do
        p(i)=p(i)-tmp*p(k)
     end do
  end do

  p(n)=p(n)/hmat(n,n) !後退代入
  do i=n-1,1,-1
     tmp=0.0d0
     do k=i+1,n
        tmp=tmp+hmat(i,k)*p(k)
     end do
     p(i)=(p(i)-tmp)/hmat(i,i)
  end do

  !---------------------------------------------------------------------------
  ! 内点計算
  allocate(pip(nip(0))) !pip(i): i番内点における複素音圧
  do i=1,nip(0) !内点のループ
     r=????? ! r=内点と点源の距離
     pip(i)=????? ! 入射波（式(89)のG(x,x^s)）、式(64)も参照のこと
     do j=1,n !要素のループ
        do ig=1,ng !ガウス積分のループ
           yxi(:)=????? !式(86)
           r=????? ! r=|x-y(ξ)| 式(87)の内点版
           dotp=???? ! dotp=(x-y(ξ))・n 式(87)の内点版
           hank=????? ! hank=H_1^(1)(k|x-y(ξ)|) 式(87)の内点版
           pip(i)=pip(i)-????? !式(89)の右辺第二項
        end do
     end do
  end do

  ! 内点の結果をplot
  ip=0
  intensity=0.d0
  open(1,file="domain.res")
  do j=0,nip(2)
     do i=0,nip(1)
        ip=ip+1
        write(1,*) xip(:,ip), real(pip(ip)), aimag(pip(ip))
     end do
     write(1,*)
  end do
  close(1)
  
  ! 結果（複素音圧の実部 = t=0における時間域の音圧分布）をreslt.pngに保存
  call system("gnuplot plot.gp")

  ! 結果（時間域の音圧）を movie.gif に保存
  call system("gnuplot mkmovie.gp")
  
end program main
