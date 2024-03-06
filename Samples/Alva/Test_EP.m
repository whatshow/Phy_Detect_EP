close all;
clear;
clc;

%%  -----------------------------------------------------------------------
%   DEBUG
%   -----------------------------------------------------------------------
%
debug = 1;                                                  % 1 do printing, 0 don't print
debug_testData = 0;

%% Param Config - Model
SNR_range = 16:2:30;                                        % SNR range
M = 16;                                                     % M-ary QAM
sym_bitnum = log2(M);                                       % Bit number in 1 M-ary modulation symbol
sympool = qammod([0: M - 1], M, "UnitAveragePower", true);  % The symbol pool to store all possible M-ary modulation symbols
sympool_real = 1/sqrt((2/3)*(M-1))*pammod([0:(2^(log2(M)/2))-1],(2^(log2(M)/2)));
sympool_var = 1;                                            % Symbol variance = Symbol Average Power                
tx_num = 12;                                                % Tx antenna number
rx_num = 12;                                                % Rx antenna number
BER_max_try_times = 1e5;                                    % The times to calculate BER to get the mean BER

%% Param Config - EP algorithm 
ep_iter_times = 10;                                         % The times to estimate the distribution of  
ep_beta = 0.9;                                              % Convex combination, usually is 0.2

%% Simulation
BERs = zeros(1, length(SNR_range));                         % BERs for every SNR
SERs = zeros(1, length(SNR_range));                         % SERs for every SNR
SERs_alva = zeros(1, length(SNR_range));                         % SERs for every SNR

for idx = 1:length(SNR_range)
    % Get current SNR
    SNR = SNR_range(idx);
    fprintf("SNR=%d\n", SNR);
    
    % Prepare the space to store all BERs during 'BER_max_try_times' times
    BER_all = zeros(1, BER_max_try_times);
    BER_all_alva = zeros(1, BER_max_try_times);
    SER_all = zeros(1, BER_max_try_times);
    SER_all_alva = zeros(1, BER_max_try_times);
    % Try several times to do average on all BERs to avoid fluctuation
    for try_times = 1:1:BER_max_try_times
        if debug_testData == 1
            nbits_len = tx_num*sym_bitnum;
            nbits = [1;0;1;0;1;0;1;1;0;0;1;0;1;0;1;1;1;0;0;0;1;1;0;0;0;0;1;1;0;0;1;0;0;0;0;1;0;0;1;0;1;1;0;1;1;0;0;1];
            H = [0.0868194352728929 + 0.113472423722110i,0.190200847087946 - 0.0888653414060696i,-0.110889806670082 - 0.362811165701612i,-0.137106721184000 + 0.392249512468358i,-0.0385385346255151 - 0.502896553677557i,0.0431292487160416 + 0.132894712220599i,-0.186035110198670 + 0.471510475902621i,0.129364166995796 - 0.326742488591624i,-0.127915488263686 + 0.460389195540822i,-0.102227284315467 + 0.210959181141221i,0.0828170703309089 + 0.168645194467220i,0.197568680602126 - 0.183797510021147i;0.147944213506750 - 0.0910870924187614i,0.240121351185112 - 0.168669947084443i,-0.0662256131762887 + 0.272422900020597i,-0.162671922994635 - 0.224095490584508i,-0.146329862701038 + 0.561534340724776i,-0.114893381106392 - 0.779818844066948i,0.125384449787720 - 0.159472165448622i,0.418154295933640 - 0.0866663304381110i,0.0697725714118238 - 0.0808872561163698i,0.224691167727634 - 0.492148533344691i,-0.264514982804913 + 0.275543496914834i,0.185241733405724 + 0.117491573263768i;-0.116154077942106 - 0.189868106270694i,-0.165910795509012 + 0.0195323657051098i,-0.102870010974272 - 0.0260610839632638i,-0.0111204998759399 + 0.360279810635718i,-0.214374002873218 - 0.402017563710757i,0.0876404712646536 - 0.307313957069288i,-0.111491767363223 + 0.166271611023925i,0.134421389850737 - 0.234296268787093i,-0.295636480092596 + 0.159243436038900i,0.201342595540077 + 0.342389208025360i,-0.335520720194665 - 0.120122859572372i,-0.118010226176369 - 0.352348366770579i;0.143239046737171 - 0.355152346213799i,-0.0199480117764521 - 0.202969056832509i,0.0185355761544619 + 0.178906323236412i,-0.0124056079995263 - 0.299222618641189i,0.0637541923891077 + 0.0437546808050478i,-0.323690231608315 - 0.283495809754432i,0.235812405975611 + 0.0199983487660457i,-0.0926992617874224 + 0.486369816659406i,0.0147609532121966 - 0.0160606041505907i,0.0390832660367674 - 0.103627927853379i,-0.229591401398610 + 0.237228093201680i,0.198956603794363 - 0.396640685721134i;-0.227402738009227 - 0.279874426534714i,0.0413692306357230 - 0.260700811459004i,0.319642837495517 - 0.215733676250825i,-0.123340111121017 - 0.0611077544195811i,-0.215548629106868 - 0.0540389775241336i,0.00339826815441289 + 0.198377429916563i,0.0413896094948915 - 0.202739828687986i,0.119559315587774 + 0.166206730055272i,-0.0268743130808401 + 0.0666154463631731i,0.304408101947947 + 0.214843070172284i,-0.0327717498854601 + 0.0192713694038527i,0.0171908889195573 - 0.0608100116640283i;0.230678010558944 + 0.0700285574303898i,-0.192360577016005 - 0.553391539787516i,-0.0387075869276703 - 0.383674649750091i,0.192446548758978 + 0.00149902882354164i,0.160275058902016 - 0.400730413087110i,0.128896027654168 - 0.274270426467850i,0.0700918495786104 - 0.0250549518185623i,-0.167909306727880 - 0.261524499503107i,0.0717459477318799 - 0.0128861071477818i,0.197751662058686 + 0.0208486889578598i,0.278440056575442 + 0.175511242992935i,-0.257864163539398 - 0.0481241419560531i;-0.203797124312630 - 0.0912920923225112i,-0.0483945532370865 - 0.0199130676260142i,-0.0456826933229819 - 0.0805227768683916i,0.157694096476678 - 0.211822636958137i,-0.0919855347600357 + 0.230730909093060i,0.0921530446837202 + 0.0901090325596421i,0.0404900412007220 - 0.194589928056918i,0.100539000646000 - 0.0511864293352531i,0.0477363862844598 - 0.101022377021375i,-0.0297713279830410 + 0.103141946674394i,0.0561650681109367 - 0.0655650542750314i,0.358371493122769 - 0.148176344556786i;0.444880596448465 + 0.0241963488186624i,0.0391524383507405 - 0.0930417724620087i,0.173087844281195 + 0.249207202775078i,0.187913350101071 - 0.0273510061286852i,-0.180306745181997 - 0.244129867130571i,-0.349622450149224 - 0.000920905014514823i,0.0225464192484733 - 0.107847778691626i,0.394721005355633 + 0.261379373688524i,0.0293975402082947 + 0.0498350164270583i,-0.127670075580915 - 0.307952008561565i,-0.122220007264528 - 0.0545179107385013i,-0.0107383691860339 - 0.303278349459960i;-0.574701889609382 - 0.166703685667401i,-0.195691381236242 - 0.231048031158418i,0.00266253698397740 - 0.149432637674194i,0.330931066048927 + 0.293293710530620i,0.171084378361112 + 0.106273610118525i,0.111070460310449 + 0.224688005705260i,0.0429727891250533 - 0.207122117005930i,-0.287824146366690 + 0.0664684972002110i,-0.0217545498798904 - 0.222707291725968i,0.228319766920741 + 0.0128436330670727i,-0.475407268996355 + 0.100989238067158i,-0.207952823587527 - 0.0232528317650572i;0.290099367192069 - 0.312884565499897i,0.0392458053410779 - 0.131409489682393i,-0.163229020141336 - 0.181043810012965i,-0.238164281785163 + 0.0791270848313183i,-0.0213267684938311 + 0.119649810490718i,-0.0888316792030517 + 0.142461745142095i,-0.149485155843498 + 0.121768644682255i,0.0306946846096202 + 0.0924039602878757i,-0.207309613047125 + 0.175305452499971i,0.0875725838173311 - 0.309396792945962i,0.308980268189903 + 0.0577306058845195i,-0.233423328048604 + 0.100607605739618i;-0.133532213852694 + 0.0151735807374821i,0.0279262949268570 - 0.00745927038162379i,0.0651571231833748 + 0.0220766383242198i,-0.130890120017566 + 0.00300674425634950i,0.0778868465217788 - 0.218234213467644i,-0.0937689216770174 + 0.193293584540364i,0.149506479597476 + 0.0508293441122704i,-0.275481886360172 - 0.143398538509297i,-0.00380290670037636 - 0.0123892080525123i,0.189269553993683 + 0.130205569861384i,0.317864070500530 - 0.108426714827986i,-0.265550549276004 + 0.376783140820528i;-0.198507270581282 + 0.0412056045094572i,0.377170185670098 - 0.480765733057030i,0.0198617223905541 + 0.201295129188095i,0.380305358834330 - 0.103561585916277i,0.153091914805322 - 0.0774740724464085i,-0.217684056138771 - 0.0595175570826921i,-0.0580713069031403 + 0.273652915292255i,0.297758716950452 + 0.227330256681059i,0.311374206647305 - 0.111792789630688i,0.0335810959474897 + 0.178081372544086i,-0.0926286161595520 - 0.144179807978595i,-0.0858698681836617 + 0.108956722269412i];
            noiseLevel = 10^(-SNR/10);
            noise = [0.0315438164023030-0.126289115870612i;-0.0328452766550056 - 0.132297806351080i;-0.115268509241838 - 0.161349443492101i;0.00524114882439885 - 0.0810088087423339i;0.173501979994306 - 0.103191131839433i;-0.0337819383892201 + 0.142212853460733i;0.0563935138759372 - 0.0148537330587628i;0.0138926432497236 + 0.0714293256621981i;0.0137469075456079 + 0.0168072546180754i;-0.0150751380809555 + 0.106910901182948i;-0.0562456431657643 + 0.189770998013110i;-0.136155704632433 - 0.0268281623924488i];
        else
            % nbits
            nbits_len = tx_num*sym_bitnum;
            nbits = randi([0 1], nbits_len, 1);
            % Channel 
            %H = 1/sqrt(2*tx_num)*randn(tx_num, rx_num) + 1/sqrt(2)*randn(tx_num, rx_num)*1j;
            H = (randn(tx_num, rx_num) + 1j*randn(tx_num, rx_num))/sqrt(2*tx_num) ;
            % Noise Creation
            noiseLevel = 10^(-SNR/10);
            noise = sqrt(noiseLevel/2) * (randn(rx_num,1) + 1j*randn(rx_num,1)) ;
            %noise = (randn(rx_num, 1) + randn(rx_num, 1)*1j)*sqrt(noiseLv/2);
        end
      
        % Create symbols
        x = qammod(nbits, M,'InputType','bit','UnitAveragePower',true);
        
        % Through AWGN channel to get y 
        y = H*x + noise;
        % EP detections
        

        y_real = [real(y);imag(y)];
        H_real = [real(H), -imag(H); imag(H), real(H)];

        [syms] = Detection_EP(sympool_real, H_real, y_real, noiseLevel/2, ep_iter_times, "Beta", ep_beta, "MinVariance", 1e-13);
        syms = [syms(1:length(syms)/2) + 1j*syms(length(syms)/2+1:end)];
        [syms_alva] = EP(y_real,H_real,noiseLevel, 'QAM', log2(M), ep_iter_times);
        % To bits
        nbits_pred = qamdemod(syms, M,'OutputType','bit','UnitAveragePower',true);
        nbits_pred_alva = qamdemod(syms_alva, M,'OutputType','bit','UnitAveragePower',true);
        % BER
        BER_all(1, try_times) = sum(nbits_pred ~= nbits)/nbits_len;  
        BER_all_alva(1, try_times) = sum(nbits_pred_alva ~= nbits)/nbits_len;  
        % SER
        %syms_pred = qammod(nbits_pred, M,'InputType','bit','UnitAveragePower',true);
        %SER_all(1, try_times) = length(find(x - syms_pred))/length(syms_pred);
        SER_all(1, try_times) = sum(x - syms > eps)/length(syms);
        SER_all_alva(1, try_times) = sum(x - syms_alva > eps)/length(syms_alva);
    end
    % do BER average
    BERs(idx) = mean(BER_all);

    SERs(idx) = mean(SER_all);
    SERs_alva(idx) = mean(SER_all_alva);
end

% plot
semilogy(SNR_range, SERs, "-ob");
hold on;
semilogy(SNR_range, SERs_alva, "-or");
grid on;
xlabel("SNR(dB)");
ylabel("SER");
ylim([10^-4, 1]);
xlim([16, 30]);
legend('EP-xin', 'EP-Alva');