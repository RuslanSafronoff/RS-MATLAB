m = 4;          % Степень расширения поля Галуа GF(2^m) - число бит в символе, m - натуральное число от 2 до 16
q = 2^m;        % Характеристика поля Галуа GF(2^m) - число различных символов поля
n = q - 1;      % Максимальная длина кодового слова Рида-Соломона - число символов в слове, n < 65536.
k = 5;          % Число информационных символов в кодовом слове
                % Расстояние Хэмминга для РС(n,k): d = n ? k + 1;
                % Расстояние Хэмминга для поризвольного (n,k): d = 2*t + 1;
                % Кратность исправляемых ошибок для РС(n,k): t = floor((n-k)/2)
J = q^k;        % Мощность РС(n,k)-кода - количество различных кодовых слов

disp(['1. Вычислим все первообразные неприводимые полиномы над GF(2^', int2str(m), '):']);
AllPrimeIrPolyDec = primpoly(m,'all'); %,'nodisplay'
AllPrimeIrPolyBin = de2bi(AllPrimeIrPolyDec, m+1, 2, 'right-msb');
disp('p(x)dec:   p0    p1    p2   ...');
disp([AllPrimeIrPolyDec, AllPrimeIrPolyBin]);

disp(' ');
disp(['2. Выберем первообразный неприводимый полином над GF(2^', int2str(m), '):']);
PrimeIrPolyDec = AllPrimeIrPolyDec(1);
PrimeIrPolyBin = de2bi(PrimeIrPolyDec, m+1, 2, 'right-msb');
disp('p(x)dec:   p0    p1    p2   ...');
disp([PrimeIrPolyDec, PrimeIrPolyBin]);

disp(' ');
disp(['3. Вычислим все элементы поля GF(2^', int2str(m), '):']);
FieldBin = gftuple([-1:q-2]', PrimeIrPolyBin);
disp('    Exp         Bin         Dec');
disp('   -----  ---------------  -----');
disp([[-1:q-2]', FieldBin, bi2de(FieldBin, 'right-msb')]);

disp(' ');
disp(['4. Вычислим коэффициенты генераторного полинома g(x) для РС(', int2str(n), ',', int2str(k), ')-кода:']);
disp('    g(X) = (X - аlpha^(b+0))(X - аlpha^(b+1))...(X - аlpha^(b+2t-1)),');
disp('    где b - целое число, как правило выбирают b = 1');
disp(['        аlpha - корень первообразного неприводимого полинома GF(2^', int2str(m), ')']);
disp('        t - кратность исправляемых ошибок, (n-k)/2');
b=1;
[GenPolyGF, t] = rsgenpoly(n,k,PrimeIrPolyDec,b);
disp(['g(X) = ', int2str(GenPolyGF.x)]);
disp(['t = ', int2str(t)]);

disp(' ');
disp(['5. Сформируем TxM - передаваемый вектор информационного слова длины k = ', int2str(k), ':']);
TxMsgWord = de2bi(floor(rand*J), k, q, 'right-msb');
TxMsgWordGF = gf(TxMsgWord, m, PrimeIrPolyDec);
disp(['TxM = ', int2str(TxMsgWord)]);

disp(' ');
disp(['6. Закодируем информационное слово TxM кодом РС(', int2str(n), ',', int2str(k), ')']);
disp(['   Сформируем TxC - вектор передаваемого кодового слова длины n = ', int2str(n), ':']);
TxCodeWordGF = rsenc(TxMsgWordGF, n, k, GenPolyGF);
disp(['TxC = ', int2str(TxCodeWordGF.x)]);

disp(' ');
disp('7. Передадим кодовое слово TxC по каналу с ошибками ...');
disp(['  Зададим E - вектор ошибок кратности t = ', int2str(t)]);
ErrWord = zeros(1, n);
for i = 1:t
    ErrPos = floor(rand*n)+1;
    ErrVal = floor(rand*q);
    ErrWord(ErrPos) = ErrVal;
end
ErrWordGF = gf(ErrWord, m, PrimeIrPolyDec);
disp(['E = ', int2str(ErrWordGF.x)]);

disp(' ');
disp('8. Вычислим RxC - принятое кодовое слово с ошибками: RxC = TxC + E');
RxCodeWordGF = TxCodeWordGF + ErrWordGF;
disp(['RxC = ', int2str(RxCodeWordGF.x)]);

disp(' ');
disp('9. Попытаемся исправить ошибки в принятом слове RxC');
disp('   Вычислим DxC и DxM  - декодированные кодовое и информационное слова');
[DxMsgWordGF, CorNumErr, DxCodeWordGF] = rsdec(RxCodeWordGF, n, k, GenPolyGF);

if CorNumErr >= 0
    disp(['Удалось исправить ', int2str(CorNumErr), ' ошибки(-ок)']);
else  
    disp('Произошла неисправимая комбинация ошибок !!!');
end
disp(['DxC = ', int2str(DxCodeWordGF.x)]);
disp(['DxM = ', int2str(DxMsgWordGF.x)]);