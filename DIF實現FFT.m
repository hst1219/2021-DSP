function X=FFT_BRO(x,log2_N)
BRO_I=BRO_Index_G(log2_N);
X=x(BRO_I+1);


function BRO_I=BRO_Index_G(log2_N)
N=2^log2_N; G=2.^[0:log2_N-1];
BRO_I=[];
for k=0:N-1
    ki=dec2bin(k,log2_N)-'0';
    BRO_I=[BRO_I ki*G.'];
end

% BF_in: 2^Lx1(L>= log2_N)
function BF_out=FFT_DIF_BF_C(BF_in,log2_N);
M=length(BF_in);
N=2^log2_N;
X=reshape(BF_in,N,M/N);
BF_no=N/2;
BF_span=N/2;
for k=1:BF_no
    A= X(k,:)+X(k+BF_span,:);
    B= X(k,:)-X(k+BF_span,:);
    X(k,:)=A;
    X(k+BF_span,:)=B;
end
BF_out=X(:);


function TW_out=FFT_DIF_TW_C(TW_in,log2_N);
M=length(TW_in);
N=2^log2_N;
X=reshape(TW_in,N,M/N);
BF_no=N/2;  BF_span=N/2;
W_N=exp(-j*2*pi/N);
for k=2:BF_no
    X(k+BF_span,:)=W_N^(k-1)*X(k+BF_span,:);
end
TW_out=X(:);



% 開始
function FFT_out=FFT_DIF_C(FFT_in,log2_N);
N=2^log2_N;
if length(FFT_in) >= N
    FFT_in=FFT_in(1:N);
else
    FFT_in=(FFT_in;zeros(N-length(FFT_in),1);
end

%一開始沒值，所以先拉到外面
X=FFT_DIF_BF_C(FFT_in,log2_N); 
X=FFT_DIF_TW_C(X,log2_N);

%中間的
for ks=log2_N-1:-1:2
    X=FFT_DIF_BF_C(X,ks);
    X=FFT_DIF_TW_C(X,ks);
end

%最後一個不用*w (有乘沒乘一樣)
X=FFT_DIF_BF_C(X,1);

% BRO
X=FFT_BRO(FFT_in,log2_N);

%結果
FFT_out=X(:);