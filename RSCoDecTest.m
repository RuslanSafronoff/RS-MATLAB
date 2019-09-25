m = 4;          % ������� ���������� ���� ����� GF(2^m) - ����� ��� � �������, m - ����������� ����� �� 2 �� 16
q = 2^m;        % �������������� ���� ����� GF(2^m) - ����� ��������� �������� ����
n = q - 1;      % ������������ ����� �������� ����� ����-�������� - ����� �������� � �����, n < 65536.
k = 5;          % ����� �������������� �������� � ������� �����
                % ���������� �������� ��� ��(n,k): d = n ? k + 1;
                % ���������� �������� ��� ������������� (n,k): d = 2*t + 1;
                % ��������� ������������ ������ ��� ��(n,k): t = floor((n-k)/2)
J = q^k;        % �������� ��(n,k)-���� - ���������� ��������� ������� ����

disp(['1. �������� ��� ������������� ������������ �������� ��� GF(2^', int2str(m), '):']);
AllPrimeIrPolyDec = primpoly(m,'all'); %,'nodisplay'
AllPrimeIrPolyBin = de2bi(AllPrimeIrPolyDec, m+1, 2, 'right-msb');
disp('p(x)dec:   p0    p1    p2   ...');
disp([AllPrimeIrPolyDec, AllPrimeIrPolyBin]);

disp(' ');
disp(['2. ������� ������������� ������������ ������� ��� GF(2^', int2str(m), '):']);
PrimeIrPolyDec = AllPrimeIrPolyDec(1);
PrimeIrPolyBin = de2bi(PrimeIrPolyDec, m+1, 2, 'right-msb');
disp('p(x)dec:   p0    p1    p2   ...');
disp([PrimeIrPolyDec, PrimeIrPolyBin]);

disp(' ');
disp(['3. �������� ��� �������� ���� GF(2^', int2str(m), '):']);
FieldBin = gftuple([-1:q-2]', PrimeIrPolyBin);
disp('    Exp         Bin         Dec');
disp('   -----  ---------------  -----');
disp([[-1:q-2]', FieldBin, bi2de(FieldBin, 'right-msb')]);

disp(' ');
disp(['4. �������� ������������ ������������� �������� g(x) ��� ��(', int2str(n), ',', int2str(k), ')-����:']);
disp('    g(X) = (X - �lpha^(b+0))(X - �lpha^(b+1))...(X - �lpha^(b+2t-1)),');
disp('    ��� b - ����� �����, ��� ������� �������� b = 1');
disp(['        �lpha - ������ �������������� ������������� �������� GF(2^', int2str(m), ')']);
disp('        t - ��������� ������������ ������, (n-k)/2');
b=1;
[GenPolyGF, t] = rsgenpoly(n,k,PrimeIrPolyDec,b);
disp(['g(X) = ', int2str(GenPolyGF.x)]);
disp(['t = ', int2str(t)]);

disp(' ');
disp(['5. ���������� TxM - ������������ ������ ��������������� ����� ����� k = ', int2str(k), ':']);
TxMsgWord = de2bi(floor(rand*J), k, q, 'right-msb');
TxMsgWordGF = gf(TxMsgWord, m, PrimeIrPolyDec);
disp(['TxM = ', int2str(TxMsgWord)]);

disp(' ');
disp(['6. ���������� �������������� ����� TxM ����� ��(', int2str(n), ',', int2str(k), ')']);
disp(['   ���������� TxC - ������ ������������� �������� ����� ����� n = ', int2str(n), ':']);
TxCodeWordGF = rsenc(TxMsgWordGF, n, k, GenPolyGF);
disp(['TxC = ', int2str(TxCodeWordGF.x)]);

disp(' ');
disp('7. ��������� ������� ����� TxC �� ������ � �������� ...');
disp(['  ������� E - ������ ������ ��������� t = ', int2str(t)]);
ErrWord = zeros(1, n);
for i = 1:t
    ErrPos = floor(rand*n)+1;
    ErrVal = floor(rand*q);
    ErrWord(ErrPos) = ErrVal;
end
ErrWordGF = gf(ErrWord, m, PrimeIrPolyDec);
disp(['E = ', int2str(ErrWordGF.x)]);

disp(' ');
disp('8. �������� RxC - �������� ������� ����� � ��������: RxC = TxC + E');
RxCodeWordGF = TxCodeWordGF + ErrWordGF;
disp(['RxC = ', int2str(RxCodeWordGF.x)]);

disp(' ');
disp('9. ���������� ��������� ������ � �������� ����� RxC');
disp('   �������� DxC � DxM  - �������������� ������� � �������������� �����');
[DxMsgWordGF, CorNumErr, DxCodeWordGF] = rsdec(RxCodeWordGF, n, k, GenPolyGF);

if CorNumErr >= 0
    disp(['������� ��������� ', int2str(CorNumErr), ' ������(-��)']);
else  
    disp('��������� ������������ ���������� ������ !!!');
end
disp(['DxC = ', int2str(DxCodeWordGF.x)]);
disp(['DxM = ', int2str(DxMsgWordGF.x)]);