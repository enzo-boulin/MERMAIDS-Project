clear all
close all

T = readtable('info_events.csv');

fs = 40.014; %Hz
T.Tag = categorical(T.Tag);
uTag = unique(T.Tag)

%%
fdata = dir('DonneesB23\*\*\*');
idx = find(arrayfun(@(x) (x.isdir==1),fdata)==0);
fdata = fdata(idx);
fdata_sig = zeros(length(fdata),1);

%%
ND = 2^(nextpow2(ceil(60*60*fs))-1);
D = nan(1,2*ND);
t = [0:2*ND-1]*1/fs;
WD = nan(1,262315);

h_ = @(x) length(strfind(T.file{k},x.name));

tagD = [];

Ddt = datetime;
kk = 0;
for k = 1:size(T,1)
    disp([k size(T,1)])
    fp = fopen(strrep(T.file{k},'data','.'),'rb');
    k_ = find(arrayfun(h_,fdata,'UniformOutput',1));
    if (~isempty(k_))
        fdata_sig(k_) = 1;
    else
        keyboard
    end
    data = fread(fp,'int32=>double')/(2^31);
    fclose(fp);
    evt_dn = datetime(T.PickedArrival(k));
    rcd_dn = datetime(T.file_start(k));
    k_evt = floor(seconds(evt_dn-rcd_dn)*fs);
    k0_evt = k_evt-ND;
    k1_evt = k_evt+ND-1;
    if (k0_evt >= 1)&(k1_evt<=length(data))
        kk = kk + 1;
        data_cn = data(k0_evt:k1_evt);
        % data_cn = (data_cn-mean(data_cn))/std(data_cn);
        D(kk,1:k1_evt-k0_evt+1) = data_cn;
        [WD(kk,:),l] = wavedec(data_cn,12,'sym8');
        Ddt(kk,1) = rcd_dn;
        Ddt(kk,2) = evt_dn;
        tagD(kk) = T.Tag(k);
    end
end
[~,idx_Ddt]=sort(Ddt(:,2),'ascend');
Ddt = Ddt(idx_Ddt,:);
D = D(idx_Ddt,:);
WD = WD(idx_Ddt,:);

%%
B = nan(1,2*ND);
WB = nan(1,262315);
Bdt = datetime;
kk = 0;
for k = 1:length(fdata)
    if (fdata_sig(k)==0)
        fp = fopen([fdata(k).folder filesep fdata(k).name],'rb');
        data = fread(fp,'int32=>double')/(2^31);
        fclose(fp);
        if (~isempty(strfind(fdata(k).name,'.')))
            Bdt_ = datetime(fdata(k).name,'InputFormat','yyyy-MM-dd''T''HH_mm_ss.SSSSSS');
        else
            Bdt_ = datetime(fdata(k).name,'InputFormat','yyyy-MM-dd''T''HH_mm_ss');
        end
        if (length(data)>2*ND)
            % data = (data-mean(data))/std(data);
            kk = kk + 1
            Bdt(kk) = Bdt_;
            B(kk,:) = data([1:2*ND]);
            [WB(kk,:),l] = wavedec(data([1:2*ND]),12,'sym8');
        end
    end
end

%%
Nh = 1024*4;
L = Nh;
f = mkf(Nh,fs);
DSPb = dsp1d(B',1/fs,L,'blackman',Nh');
DSPd = dsp1d(D',1/fs,L,'blackman',Nh');

DSPb_Ddt = interp1(Bdt,DSPb,Ddt(:,2),'lin');

pDSPd = DSPd(:,Nh/2+1:end);
pDSPb_Ddt = DSPb_Ddt(:,Nh/2+1:end);


%%
H = max(0,mean(DSPd,1)-mean(DSPb_Ddt,1,'omitmissing'))./(mean(DSPd,1));
H(abs(f)>0.1) = 0;

idx_f = find(abs(f)>2);
h = real(ftm2ri(H));
Dh = nan(size(D));
Dh_ = nan(size(D));
H_ = [];
for k = 1:size(D,1)
    Dh(k,:) = conv(D(k,:),h,'same');
    H_(k,:) = max(0,DSPd(k,:)-nanmean(DSPb_Ddt))./(DSPd(k,:));
    H_(k,idx_f) = 0;
    h_ = real(ftm2ri(H_(k,:)));
    Dh_(k,:) = conv(D(k,:),h_,'same');
end
DSPdh = dsp1d(Dh',1/fs,L,'blackman',Nh)';

Bh = nan(size(B));
for k = 1:size(B,1)
    Bh(k,:) = conv(B(k,:),h,'same');
end
DSPbh = dsp1d(Bh',1/fs,L,'blackman',Nh)';

%%
levels = 1:12;
longs = l;
first = cumsum(longs)+1;
first = first(end-2:-1:1);
longs = longs(end-1:-1:2);
last  = first+longs-1;
nblev = length(levels);
idx   = cell(1,nblev);
varWB = nan(size(B,1),nblev);
varWD = nan(size(D,1),nblev);
varWB_Ddt = interp1(Bdt,varWB,Ddt(:,2),'lin');

WD_ = 0*WD;

for j = 1:nblev
    k = levels(j);
    idx{j} = (first(k):last(k));
    varWB(:,j) = var(WB(:,idx{j}),0,2,'omitmissing');
    varWD(:,j) = var(WD(:,idx{j}),0,2,'omitmissing');
    varWB_Ddt_j = interp1(Bdt,varWB(:,j),Ddt(:,2),'lin');
    WD_j = WD(:,idx{j});
    % WD_j = WD_j.*max(0,stdWD(:,j)-stdWB_Ddt_j)./stdWD(:,j);
    % WD_j = WD_j.*max(0,varWD(:,j)-varWB_Ddt_j)./varWD(:,j);
    WD_j(abs(WD_j)<1.5*varWB_Ddt_j.^0.5) = 0;
    WD_(:,idx{j}) = WD_j;
end

D_ = zeros(size(D));
for k = 1:size(D,1)
    D_(k,:) = waverec(WD_(k,:),l,'sym8');
end
