%% function ESN_population_coding_all_saccades
function ESN_population_coding_all_saccades
%% build and save ALL_PCELL_COMPRESSED_DATA
path_data_monkey_sorted = uigetdir;

% build pCell_list, this is a hard coded cell with the id of all of the pCells and the bundles
pCell_list_1 = build_pCell_list_Mirza_pre201906();
pCell_list_2 = build_pCell_list_Mirza_post201906();
pCell_list_3 = build_pCell_list_Ramon();
pCell_list = vertcat(pCell_list_1, pCell_list_2, pCell_list_3);
%%
% build ALL_PCELL_COMPRESSED_DATA, open plot data and put them together
ALL_PCELL_all_saccades = build_ALL_PCELL_all_saccades(pCell_list, path_data_monkey_sorted);

num_pCells = length(pCell_list);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~exist([path_data_monkey_sorted filesep ALL_PCELL_name], 'dir')
    mkdir([path_data_monkey_sorted filesep ALL_PCELL_name]);
end

% save ALL_PCELL_COMPRESSED_DATA
save([path_data_monkey_sorted filesep ALL_PCELL_name filesep 'ALL_PCELL_all_saccades.mat'],...
    'ALL_PCELL_all_saccades',...
    '-v7.3')

%%
clearvars -except ALL_PCELL_all_saccades
ALL_PCELL_all_saccades_tuned = build_ALL_PCELL_all_saccades_tuned(ALL_PCELL_all_saccades);
ALL_PCELL_all_saccades_tuned_amp = avg_over_angle(ALL_PCELL_all_saccades_tuned);
ALL_PCELL_all_saccades_tuned_ang = avg_over_amplitude(ALL_PCELL_all_saccades_tuned);
%
ALL_PCELL_all_saccades_tuned_ang_amp = avg_over_angle(ALL_PCELL_all_saccades_tuned_ang);
ALL_PCELL_all_saccades_tuned_amp_ang = avg_over_amplitude(ALL_PCELL_all_saccades_tuned_amp);

%



end

%% function build_pCell_list_Mirza_pre201906
function pCell_list_Mirza_pre201906 = build_pCell_list_Mirza_pre201906()
%% build pCell_list
pCell_list_Mirza_pre201906 = cell(0,10);

% pCell_list(0) = {'', '', '', '', '', '', '', '', '', ''};

pCell_list_Mirza_pre201906( 1, 1:2) = {'190326_180641_03_sorted_RSh', '190326_182755_03_sorted_RSh'};
pCell_list_Mirza_pre201906( 2, 1:2) = {'190326_180641_04_sorted_RSh', '190326_182755_04_sorted_RSh'};
pCell_list_Mirza_pre201906( 3, 1:3) = {'190327_135803_05_sorted_RSh', '190327_143005_05_sorted_RSh', '190327_150203_05_sorted_RSh'};
pCell_list_Mirza_pre201906( 4, 1:3) = {'190327_135803_06_sorted_RSh', '190327_143005_06_sorted_RSh', '190327_150203_06_sorted_RSh'};
pCell_list_Mirza_pre201906( 5, 1:2) = {'190327_152414_03_sorted_RSh', '190327_153410_03_sorted_RSh'};
pCell_list_Mirza_pre201906( 6, 1:1) = {'190329_133538_05_sorted_RSh'};
pCell_list_Mirza_pre201906( 7, 1:1) = {'190329_133538_03_sorted_RSh'};
pCell_list_Mirza_pre201906( 8, 1:1) = {'190329_150427_03_sorted_RSh'};
pCell_list_Mirza_pre201906( 9, 1:1) = {'190329_150427_07_sorted_RSh'};
pCell_list_Mirza_pre201906(10, 1:1) = {'190402_140923_05_sorted_RSh'};

pCell_list_Mirza_pre201906(11, 1:1) = {'190402_150756_05_sorted_RSh'};
pCell_list_Mirza_pre201906(12, 1:2) = {'190403_141313_04_sorted_RSh', '190403_143842_04_sorted_RSh'};
pCell_list_Mirza_pre201906(13, 1:3) = {'190403_151506_03_sorted_RSh', '190403_153306_03_sorted_RSh', '190403_155643_03_sorted_RSh'};
pCell_list_Mirza_pre201906(14, 1:1) = {'190404_110707_07_sorted_RSh'};
pCell_list_Mirza_pre201906(15, 1:3) = {'190408_130716_04_sorted_ESN', '190408_131546_06_sorted_ESN', '190408_133026_06_sorted_ESN'};
pCell_list_Mirza_pre201906(16, 1:1) = {'190409_143035_04_sorted_ESN'};
pCell_list_Mirza_pre201906(17, 1:1) = {'190409_150223_03_sorted_ESN'};
pCell_list_Mirza_pre201906(18, 1:1) = {'190409_150223_07_sorted_ESN'};
pCell_list_Mirza_pre201906(19, 1:2) = {'190409_151632_03_sorted_ESN', '190409_152524_03_sorted_ESN'};
pCell_list_Mirza_pre201906(20, 1:3) = {'190409_152524_07_sorted_ESN', '190409_154839_07_sorted_ESN', '190409_155447_07_sorted_ESN'};

pCell_list_Mirza_pre201906(21, 1:3) = {'190410_140235_04_sorted_RSh', '190410_144541_04_sorted_RSh', '190410_150708_04_sorted_RSh'};
pCell_list_Mirza_pre201906(22, 1:1) = {'190411_124025_04_sorted_RSh'};
pCell_list_Mirza_pre201906(23, 1:2) = {'190411_130256_03_sorted_RSh', '190411_133808_03_sorted_ESN'};
pCell_list_Mirza_pre201906(24, 1:1) = {'190411_130256_07_sorted_ESN'};
pCell_list_Mirza_pre201906(25, 1:5) = {'190412_115230_05_sorted_RSh', '190412_120702_05_sorted_RSh', '190412_122025_05_sorted_RSh', '190412_123356_05_sorted_RSh', '190412_125020_05_sorted_RSh'};
pCell_list_Mirza_pre201906(26, 1:6) = {'190412_115230_01_sorted_RSh', '190412_120702_03_sorted_RSh', '190412_122025_03_sorted_RSh', '190412_123356_03_sorted_RSh', '190412_125020_01_sorted_RSh', '190412_131235_01_sorted_RSh'};
pCell_list_Mirza_pre201906(27, 1:4) = {'190415_135953_04_sorted_RSh', '190415_142550_05_sorted_ESN', '190415_144515_05_sorted_RSh', '190415_145524_05_sorted_RSh'};
pCell_list_Mirza_pre201906(28, 1:2) = {'190415_153818_07_sorted_RSh', '190415_160746_07_sorted_RSh'};
pCell_list_Mirza_pre201906(29, 1:2) = {'190416_142432_02_sorted_RSh', '190416_144104_02_sorted_RSh'};
pCell_list_Mirza_pre201906(30, 1:2) = {'190422_132817_04_sorted_RSh', '190422_133419_04_sorted_RSh'};

pCell_list_Mirza_pre201906(31, 1:2) = {'190422_143306_04_sorted_RSh', '190422_144948_05_sorted_RSh'};
pCell_list_Mirza_pre201906(32, 1:1) = {'190422_144948_07_sorted_RSh'};
pCell_list_Mirza_pre201906(33, 1:3) = {'190423_132324_01_sorted_ESN', '190423_134600_01_sorted_ESN', '190423_142023_01_sorted_ESN'};
pCell_list_Mirza_pre201906(34, 1:3) = {'190424_151111_04_sorted_RSh', '190424_153357_05_sorted_RSh', '190424_160407_05_sorted_RSh'};
pCell_list_Mirza_pre201906(35, 1:2) = {'190426_120355_05_sorted_RSh', '190426_122359_05_sorted_RSh'};
pCell_list_Mirza_pre201906(36, 1:2) = {'190426_120355_06_sorted_RSh', '190426_122359_06_sorted_RSh'};
pCell_list_Mirza_pre201906(37, 1:3) = {'190426_125234_07_sorted_RSh', '190426_131912_07_sorted_RSh', '190426_134533_07_sorted_RSh'};
pCell_list_Mirza_pre201906(38, 1:3) = {'190426_125234_01_sorted_RSh', '190426_131912_01_sorted_RSh', '190426_134533_01_sorted_RSh'};
pCell_list_Mirza_pre201906(39, 1:3) = {'190429_133647_04_sorted_RSh', '190429_135343_05_sorted_RSh', '190429_141529_05_sorted_RSh'};
pCell_list_Mirza_pre201906(40, 1:2) = {'190429_135343_06_sorted_RSh', '190429_141529_06_sorted_RSh'};

pCell_list_Mirza_pre201906(41, 1:3) = {'190429_142903_02_sorted_RSh', '190429_145043_02_sorted_JSP', '190429_150541_02_sorted_JSP'};
pCell_list_Mirza_pre201906(42, 1:3) = {'190429_142903_03_sorted_RSh', '190429_145043_03_sorted_RSh', '190429_150541_03_sorted_RSh'};
pCell_list_Mirza_pre201906(43, 1:3) = {'190429_142903_07_sorted_RSh', '190429_145043_07_sorted_RSh', '190429_150541_07_sorted_RSh'};
pCell_list_Mirza_pre201906(44, 1:2) = {'190430_132054_05_sorted_ESN', '190430_132215_05_sorted_ESN'};
pCell_list_Mirza_pre201906(45, 1:2) = {'190430_134129_03_sorted_RSh', '190430_140812_03_sorted_RSh'};
pCell_list_Mirza_pre201906(46, 1:1) = {'190501_141616_03_sorted_RSh'};
pCell_list_Mirza_pre201906(47, 1:1) = {'190515_135233_04_sorted_ESN'};
pCell_list_Mirza_pre201906(48, 1:4) = {'190516_144904_02_sorted_RSh', '190516_150909_02_sorted_RSh', '190516_153805_02_sorted_RSh', '190516_155125_02_sorted_RSh'};
pCell_list_Mirza_pre201906(49, 1:2) = {'190516_144904_07_sorted_RSh', '190516_150909_07_sorted_RSh'};
pCell_list_Mirza_pre201906(50, 1:1) = {'190517_113915_05_sorted_RSh'};

pCell_list_Mirza_pre201906(51, 1:2) = {'190517_121917_05_sorted_RSh', '190517_122246_05_sorted_RSh'};

end

%% function build_pCell_list_Mirza_post201906
function pCell_list_Mirza_post201906 = build_pCell_list_Mirza_post201906()
%% build pCell_list
pCell_list_Mirza_post201906 = cell(0,10);

% pCell_list(0) = {'', '', '', '', '', '', '', '', '', ''};
pCell_list_Mirza_post201906( 1, 1:5) = {'190724_141024_04_sorted_PGH', '190724_142939_05_sorted_PGH', '190724_145556_06_sorted_PGH', '190724_151910_06_sorted_PGH', '190724_154019_06_sorted_PGH'};
pCell_list_Mirza_post201906( 2, 1:4) = {'190730_143843_01_sorted_PGH', '190730_150548_01_sorted_PGH', '190730_154007_01_sorted_PGH', '190730_161519_01_sorted_PGH'};
pCell_list_Mirza_post201906( 3, 1:1) = {'190801_144133_05_sorted_RSh'};
pCell_list_Mirza_post201906( 4, 1:1) = {'190801_151921_04_sorted_RSh'};
pCell_list_Mirza_post201906( 5, 1:1) = {'190801_161911_05_sorted_RSh'};
pCell_list_Mirza_post201906( 6, 1:2) = {'190812_140139_03_sorted_RSh', '190812_142255_03_sorted_RSh'};
pCell_list_Mirza_post201906( 7, 1:3) = {'190812_153354_02_sorted_RSh', '190812_160450_02_sorted_RSh', '190812_163649_02_sorted_RSh'};
pCell_list_Mirza_post201906( 8, 1:3) = {'190815_161258_08_sorted_PGH', '190815_165738_08_sorted_PGH', '190815_171400_08_sorted_PGH'};
pCell_list_Mirza_post201906( 9, 1:2) = {'190918_115126_04_sorted_RSh', '190918_122153_04_sorted_RSh'};
pCell_list_Mirza_post201906(10, 1:2) = {'190919_095155_02_sorted_RSh', '190919_100229_02_sorted_RSh'};

pCell_list_Mirza_post201906(11, 1:3) = {'190923_143725_02_sorted_PGH', '190923_150658_02_sorted_PGH', '190923_154446_02_sorted_PGH'};
pCell_list_Mirza_post201906(12, 1:3) = {'190925_145115_02_sorted_PGH', '190925_152228_02_sorted_PGH', '190925_155449_02_sorted_PGH'};
pCell_list_Mirza_post201906(13, 1:1) = {'190927_135905_01_sorted_RSh'};
pCell_list_Mirza_post201906(14, 1:1) = {'190927_135905_03_sorted_RSh'};

end

%% function build_pCell_list_Ramon
function pCell_list_Ramon = build_pCell_list_Ramon()
%% build pCell_list
pCell_list_Ramon = cell(0,10);

% pCell_list(0) = {'', '', '', '', '', '', '', '', '', ''};
pCell_list_Ramon( 1, 1:3) = {'190829_131625_04_sorted_ESN', '190829_132447_04_sorted_ESN', '190829_133438_04_sorted_ESN'};
pCell_list_Ramon( 2, 1:2) = {'190830_114101_02_sorted_ESN', '190830_114850_04_sorted_ESN'};
pCell_list_Ramon( 3, 1:2) = {'190830_114101_04_sorted_ESN', '190830_114850_02_sorted_ESN'};
pCell_list_Ramon( 4, 1:2) = {'190903_160127_04_sorted_ESN', '190903_163920_04_sorted_ESN'};
pCell_list_Ramon( 5, 1:2) = {'190903_160127_01_sorted_ESN', '190903_163920_01_sorted_ESN'};
pCell_list_Ramon( 6, 1:1) = {'190905_154933_03_sorted_ESN'};
pCell_list_Ramon( 7, 1:1) = {'190906_125220_01_sorted_ESN'};
pCell_list_Ramon( 8, 1:2) = {'190917_144329_04_sorted_ESN', '190917_155509_04_sorted_ESN'};
pCell_list_Ramon( 9, 1:2) = {'190917_150252_02_sorted_ESN', '190917_152925_02_sorted_ESN'};
pCell_list_Ramon(10, 1:2) = {'190918_152830_02_sorted_ESN', '190918_154008_02_sorted_ESN'};

pCell_list_Ramon(11, 1:1) = {'190919_154206_02_sorted_ESN'};
pCell_list_Ramon(12, 1:2) = {'190925_113310_02_sorted_ESN', '190925_121201_02_sorted_ESN'};
pCell_list_Ramon(13, 1:1) = {'191010_122539_05_sorted_ESN'};
pCell_list_Ramon(14, 1:1) = {'191010_125057_02_sorted_ESN'};
pCell_list_Ramon(15, 1:3) = {'191101_130349_03_sorted_ESN', '191101_133125_03_sorted_ESN', '191101_140646_03_sorted_ESN'};
pCell_list_Ramon(16, 1:3) = {'191101_130349_04_sorted_ESN', '191101_133125_04_sorted_ESN', '191101_140646_04_sorted_ESN'};
pCell_list_Ramon(17, 1:3) = {'191101_130349_05_sorted_ESN', '191101_133125_05_sorted_ESN', '191101_140646_05_sorted_ESN'};
pCell_list_Ramon(18, 1:2) = {'191104_123450_02_sorted_ESN', '191104_125346_02_sorted_ESN'};
pCell_list_Ramon(19, 1:1) = {'191104_131449_03_sorted_ESN'};
pCell_list_Ramon(20, 1:2) = {'191108_113612_04_sorted_ESN', '191108_114524_04_sorted_ESN'};

pCell_list_Ramon(21, 1:1) = {'191108_124201_03_sorted_ESN'};
pCell_list_Ramon(22, 1:1) = {'191108_124201_04_sorted_ESN'};
pCell_list_Ramon(23, 1:1) = {'191108_135053_07_sorted_ESN'};
pCell_list_Ramon(24, 1:4) = {'191111_131030_04_sorted_ESN', '191111_132846_04_sorted_ESN', '191111_135712_04_sorted_ESN', '191111_142741_04_sorted_ESN'};
pCell_list_Ramon(25, 1:1) = {'191111_131030_04_2_sorted_ESN'};
pCell_list_Ramon(26, 1:8) = {'191115_115149_04_sorted_RSh', '191115_120043_04_sorted_RSh', '191115_120257_04_sorted_RSh', '191115_122302_04_sorted_RSh', '191115_122724_04_sorted_RSh', '191115_124011_04_sorted_RSh', '191115_125241_04_sorted_RSh', '191115_125423_04_sorted_RSh'};
pCell_list_Ramon(27, 1:3) = {'191118_143441_05_sorted_RSh', '191118_145837_05_sorted_RSh', '191118_154644_05_sorted_RSh'};
pCell_list_Ramon(28, 1:3) = {'191118_145837_01_sorted_RSh', '191118_152411_01_sorted_RSh', '191118_154644_01_sorted_RSh'};
pCell_list_Ramon(29, 1:3) = {'191120_141613_02_sorted_RSh', '191120_144254_02_sorted_RSh', '191120_145548_02_sorted_RSh'};
pCell_list_Ramon(30, 1:1) = {'191204_115824_03_sorted_RSh'};

pCell_list_Ramon(31, 1:1) = {'191204_131009_01_sorted_RSh'};
pCell_list_Ramon(32, 1:2) = {'200109_110904_04_sorted_RSh', '200109_112601_04_sorted_RSh'};
pCell_list_Ramon(33, 1:1) = {'200109_125642_03_sorted_RSh'};
pCell_list_Ramon(34, 1:3) = {'200122_125828_07_sorted_RSh', '200122_133020_07_sorted_RSh', '200122_135404_07_sorted_RSh'};
pCell_list_Ramon(35, 1:2) = {'200122_133020_04_sorted_RSh', '200122_135404_04_sorted_RSh'};
pCell_list_Ramon(36, 1:4) = {'200123_124808_07_sorted_PGH', '200123_131422_07_sorted_PGH', '200123_133835_07_sorted_PGH', '200123_140924_07_sorted_PGH'};
pCell_list_Ramon(37, 1:4) = {'200128_151251_04_sorted_RSh', '200128_151937_06_sorted_RSh', '200128_153045_06_sorted_RSh', '200128_154444_01_sorted_RSh'};
pCell_list_Ramon(38, 1:1) = {'200128_163856_04_sorted_RSh'};
pCell_list_Ramon(39, 1:4) = {'200129_094755_04_sorted_PGH', '200129_124934_04_sorted_PGH', '200129_131018_04_sorted_PGH', '200129_133230_04_sorted_PGH'};
pCell_list_Ramon(40, 1:3) = {'200204_162413_04_sorted_RSh', '200204_163427_04_sorted_RSh', '200204_164913_03_sorted_RSh'};

pCell_list_Ramon(41, 1:1) = {'200211_162022_03_sorted_RSh'};
pCell_list_Ramon(42, 1:5) = {'200303_153615_03_sorted_RSh', '200303_155436_03_sorted_RSh', '200303_163555_03_sorted_RSh', '200303_165850_03_sorted_RSh', '200303_172913_03_sorted_RSh'};
pCell_list_Ramon(43, 1:3) = {'200814_113452_05_sorted_RSh', '200814_114000_05_sorted_RSh', '200814_114425_05_sorted_RSh'};
pCell_list_Ramon(44, 1:3) = {'200914_122203_04_sorted_RSh', '200914_124954_04_sorted_RSh', '200914_132356_04_sorted_RSh'};
pCell_list_Ramon(45, 1:3) = {'200916_130045_02_sorted_PGH', '200916_133908_02_sorted_PGH', '200916_142429_02_sorted_PGH'};

end

%% function build_ALL_PCELL_all_saccades
function ALL_PCELL_all_saccades = build_ALL_PCELL_all_saccades(pCell_list, path_data_monkey_sorted)
%% init vars
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
clearvars ALL_PCELL_COMPRESSED_DATA
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    clearvars TRAIN_DATA_recordings
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name = pCell_list{counter_pCell, counter_recording};
        year_ = file_name(1:2);
        month_ = file_name(3:4);
        day_ = file_name(5:6);
        hour_ = file_name(8:9);
        minute_ = file_name(10:11);
        second_ = file_name(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_data = ['analyzed_data' filesep];
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_data];
        %% load plot_data
        TRAIN_DATA_recording = build_TRAIN_DATA(file_path, file_name);
        TRAIN_DATA_recording.file_name = file_name;
        TRAIN_DATA_recording.file_path = file_path;
        if TRAIN_DATA_recording.file_name(18) == 's'
            TRAIN_DATA_recording.id          = TRAIN_DATA_recording.file_name(1:16);
        elseif TRAIN_DATA_recording.file_name(18) == '2'
            TRAIN_DATA_recording.id          = TRAIN_DATA_recording.file_name(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        TRAIN_DATA_recordings(counter_recording) = TRAIN_DATA_recording;
    end
    %% 
    TRAIN_DATA_cell = struct;
    if num_recording > 1
        
        for counter_variable  = 1 : length(variable_list)
        for counter_indType   = 1 : length(indType_list)
        for counter_spikeType = 1 : length(spikeType_list)
            variable_name = variable_list{counter_variable};
            indType_name = indType_list{counter_indType};
            spikeType_name = spikeType_list{counter_spikeType};
            num_amp_bin = size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_train_' indType_name]),1);
            num_ang_bin = size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_train_' indType_name]),2);

            for counter_amp = 1 : num_amp_bin
            for counter_ang = 1 : num_ang_bin
                train_data_recordings = zeros(size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang}));
                velocity_data_recordings = zeros(size(TRAIN_DATA_recordings(1).(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang}));
                num_saccades_recordings = 0;
                for counter_recording = 1 : 1 : num_recording
                    train_data = TRAIN_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
                    velocity_data = TRAIN_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
                    num_saccades = TRAIN_DATA_recordings(counter_recording).(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
                    train_data_recordings = train_data_recordings + (train_data * num_saccades);
                    velocity_data_recordings = velocity_data_recordings + (velocity_data * num_saccades);
                    num_saccades_recordings = num_saccades_recordings + num_saccades;
                end
                if num_saccades_recordings == 0
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades_recordings;
                else
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data_recordings ./ num_saccades_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data_recordings ./ num_saccades_recordings;
                    TRAIN_DATA_cell.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades_recordings;
                end
                
            end
            end
        end
        end
        end
        % concatenate cell id
        TRAIN_DATA_cell.id = cell(num_recording, 1);
        TRAIN_DATA_cell.file_name = cell(num_recording, 1);
        TRAIN_DATA_cell.file_path = cell(num_recording, 1);
        for counter_recording = 1 : 1 : num_recording
            TRAIN_DATA_cell.id{counter_recording, 1} = TRAIN_DATA_recordings(counter_recording).id;
            TRAIN_DATA_cell.file_name{counter_recording, 1} = TRAIN_DATA_recordings(counter_recording).file_name;
            TRAIN_DATA_cell.file_path{counter_recording, 1} = TRAIN_DATA_recordings(counter_recording).file_path;
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % concatenate Neural_Properties
%         if ~isfield(TRAIN_DATA_cell, 'Neural_Properties')
%             TRAIN_DATA_cell.Neural_Properties = struct();
%         end
%         if ~isfield(TRAIN_DATA_cell.Neural_Properties, 'SS_num')
%             TRAIN_DATA_cell.Neural_Properties.SS_num = [];
%             TRAIN_DATA_cell.Neural_Properties.SS_duration = [];
%             TRAIN_DATA_cell.Neural_Properties.SS_firing_rate = [];
%             TRAIN_DATA_cell.Neural_Properties.SS_time = [];
%             TRAIN_DATA_cell.Neural_Properties.SS_waveform = [];
%         end
%         if ~isempty(plot_data_raw.Neural_Properties.CH_sorted.SS_data.SS_time)
%             TRAIN_DATA_cell.Neural_Properties.SS_num(counter_recording, 1) = ...
%                 length(plot_data_raw.Neural_Properties.CH_sorted.SS_data.SS_time);
%             TRAIN_DATA_cell.Neural_Properties.SS_duration(counter_recording, 1) = ...
%                 (plot_data_raw.Neural_Properties.CH_sorted.SS_data.SS_time(end)) - (plot_data_raw.Neural_Properties.CH_sorted.SS_data.SS_time(1));
%             TRAIN_DATA_cell.Neural_Properties.SS_firing_rate(counter_recording, 1) = ...
%                 TRAIN_DATA_cell.Neural_Properties.SS_num(counter_recording, 1) / TRAIN_DATA_cell.Neural_Properties.SS_duration(counter_recording, 1);
%             SS_time_plot_data_raw = plot_data_raw.Neural_Properties.CH_sorted.SS_data.SS_time;
%             SS_time_TRAIN_DATA_cell = TRAIN_DATA_cell.Neural_Properties.SS_time;
%             SS_time_TRAIN_DATA_cell = vertcat(SS_time_TRAIN_DATA_cell, SS_time_plot_data_raw);
%             TRAIN_DATA_cell.Neural_Properties.SS_time = SS_time_TRAIN_DATA_cell;
%             SS_waveform_plot_data_raw = plot_data_raw.Neural_Properties.CH_sorted.SS_data.SS_waveform;
%             SS_waveform_TRAIN_DATA_cell = TRAIN_DATA_cell.Neural_Properties.SS_waveform;
%             SS_waveform_TRAIN_DATA_cell = vertcat(SS_waveform_TRAIN_DATA_cell, SS_waveform_plot_data_raw);
%             TRAIN_DATA_cell.Neural_Properties.SS_waveform = SS_waveform_TRAIN_DATA_cell;
%         else
%             TRAIN_DATA_cell.Neural_Properties.SS_num(counter_recording, 1) = 0;
%             TRAIN_DATA_cell.Neural_Properties.SS_duration(counter_recording, 1) = ...
%                 (plot_data_raw.Neural_Properties.CH_sorted.CS_data.CS_time(end)) - (plot_data_raw.Neural_Properties.CH_sorted.CS_data.CS_time(1));
%             TRAIN_DATA_cell.Neural_Properties.SS_firing_rate(counter_recording, 1) = 0;
%         end
%         
%         if ~isfield(TRAIN_DATA_cell.Neural_Properties, 'CS_num')
%             TRAIN_DATA_cell.Neural_Properties.CS_num = [];
%             TRAIN_DATA_cell.Neural_Properties.CS_firing_rate = [];
%             TRAIN_DATA_cell.Neural_Properties.CS_time = [];
%             TRAIN_DATA_cell.Neural_Properties.CS_waveform = [];
%         end
%         if ~isempty(plot_data_raw.Neural_Properties.CH_sorted.CS_data.CS_time)
%             TRAIN_DATA_cell.Neural_Properties.CS_num(counter_recording, 1) = ...
%                 length(plot_data_raw.Neural_Properties.CH_sorted.CS_data.CS_time);
%             TRAIN_DATA_cell.Neural_Properties.CS_firing_rate(counter_recording, 1) = ...
%                 TRAIN_DATA_cell.Neural_Properties.CS_num(counter_recording, 1) / TRAIN_DATA_cell.Neural_Properties.SS_duration(counter_recording, 1);
%             CS_time_plot_data_raw = plot_data_raw.Neural_Properties.CH_sorted.CS_data.CS_time;
%             CS_time_TRAIN_DATA_cell = TRAIN_DATA_cell.Neural_Properties.CS_time;
%             CS_time_TRAIN_DATA_cell = vertcat(CS_time_TRAIN_DATA_cell, CS_time_plot_data_raw);
%             TRAIN_DATA_cell.Neural_Properties.CS_time = CS_time_TRAIN_DATA_cell;
%             CS_waveform_plot_data_raw = plot_data_raw.Neural_Properties.CH_sorted.CS_data.CS_waveform;
%             CS_waveform_TRAIN_DATA_cell = TRAIN_DATA_cell.Neural_Properties.CS_waveform;
%             CS_waveform_TRAIN_DATA_cell = vertcat(CS_waveform_TRAIN_DATA_cell, CS_waveform_plot_data_raw);
%             TRAIN_DATA_cell.Neural_Properties.CS_waveform = CS_waveform_TRAIN_DATA_cell;
%         else
%             TRAIN_DATA_cell.Neural_Properties.CS_num(counter_recording, 1) = 0;
%             TRAIN_DATA_cell.Neural_Properties.CS_firing_rate(counter_recording, 1) = 0;
%         end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        TRAIN_DATA_cell.Neural_Properties.SS_num = 0;
        TRAIN_DATA_cell.Neural_Properties.SS_duration = 0;
        TRAIN_DATA_cell.Neural_Properties.SS_firing_rate = 0;
        TRAIN_DATA_cell.Neural_Properties.SS_time = [];
        TRAIN_DATA_cell.Neural_Properties.SS_waveform = [];
        TRAIN_DATA_cell.Neural_Properties.CS_num = 0;
        TRAIN_DATA_cell.Neural_Properties.CS_firing_rate = 0;
        TRAIN_DATA_cell.Neural_Properties.CS_time = [];
        TRAIN_DATA_cell.Neural_Properties.CS_waveform = [];
        for counter_recording = 1 : 1 : num_recording
            TRAIN_DATA_cell.Neural_Properties.SS_num = TRAIN_DATA_cell.Neural_Properties.SS_num + ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_num;
            TRAIN_DATA_cell.Neural_Properties.SS_duration = TRAIN_DATA_cell.Neural_Properties.SS_duration + ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_duration;
            TRAIN_DATA_cell.Neural_Properties.SS_time = vertcat(TRAIN_DATA_cell.Neural_Properties.SS_time ,...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_time);
            TRAIN_DATA_cell.Neural_Properties.SS_waveform = vertcat(TRAIN_DATA_cell.Neural_Properties.SS_waveform, ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.SS_waveform);

            TRAIN_DATA_cell.Neural_Properties.CS_num = TRAIN_DATA_cell.Neural_Properties.CS_num + ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.CS_num;
            TRAIN_DATA_cell.Neural_Properties.CS_time = vertcat(TRAIN_DATA_cell.Neural_Properties.CS_time, ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.CS_time);
            TRAIN_DATA_cell.Neural_Properties.CS_waveform = vertcat(TRAIN_DATA_cell.Neural_Properties.CS_waveform, ...
                TRAIN_DATA_recordings(counter_recording).Neural_Properties.CS_waveform);
            % concatenate Neural_Properties.Corr_data
            if ~isfield(TRAIN_DATA_cell.Neural_Properties, 'Corr_data')
                TRAIN_DATA_cell.Neural_Properties.Corr_data = struct();
            end
            field_names_Corr_data = fieldnames(TRAIN_DATA_recordings(counter_recording).Neural_Properties.Corr_data);
            for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
                field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
                data_field_name_Corr_data = TRAIN_DATA_recordings(counter_recording).Neural_Properties.Corr_data.(field_name_Corr_data);
                if ~isfield(TRAIN_DATA_cell.Neural_Properties.Corr_data, field_name_Corr_data)
                    TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = [];
                end
                data_field_name_cell = TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data);
                data_field_name_cell = vertcat(data_field_name_cell, data_field_name_Corr_data);
                TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = data_field_name_cell;
            end
        end
        TRAIN_DATA_cell.Neural_Properties.SS_firing_rate = ...
                TRAIN_DATA_cell.Neural_Properties.SS_num / TRAIN_DATA_cell.Neural_Properties.SS_duration;
        TRAIN_DATA_cell.Neural_Properties.CS_firing_rate = ...
                TRAIN_DATA_cell.Neural_Properties.CS_num / TRAIN_DATA_cell.Neural_Properties.SS_duration;
        TRAIN_DATA_cell.Neural_Properties.SS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.SS_waveform, 1);
        TRAIN_DATA_cell.Neural_Properties.CS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.CS_waveform, 1);
        field_names_Corr_data = fieldnames(TRAIN_DATA_cell.Neural_Properties.Corr_data);
        for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
            field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
            data_field_name_Corr_data = TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data);
            TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = nanmean(data_field_name_Corr_data, 1);
        end
    else
        TRAIN_DATA_cell = TRAIN_DATA_recordings;
        TRAIN_DATA_cell.Neural_Properties.SS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.SS_waveform, 1);
        TRAIN_DATA_cell.Neural_Properties.CS_waveform = nanmean(TRAIN_DATA_cell.Neural_Properties.CS_waveform, 1);
        field_names_Corr_data = fieldnames(TRAIN_DATA_cell.Neural_Properties.Corr_data);
        for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
            field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
            data_field_name_Corr_data = TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data);
            TRAIN_DATA_cell.Neural_Properties.Corr_data.(field_name_Corr_data) = nanmean(data_field_name_Corr_data, 1);
        end
    end
    
    %% Save TRAIN_DATA_cell into ALL_PCELL_COMPRESSED_DATA
    ALL_PCELL_COMPRESSED_DATA(counter_pCell) = TRAIN_DATA_cell;
    fprintf(' --> Completed. \n')
end

%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades', ' ... ']);
clearvars ALL_PCELL_all_saccades
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_COMPRESSED_DATA(1).(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_COMPRESSED_DATA(1).(variable_name).([spikeType_name '_train_' indType_name]),2);

    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        for counter_pCell = 1 : 1 : num_pCells
            train_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_COMPRESSED_DATA(counter_pCell).(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang);
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = train_data;
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = velocity_data;
            ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = num_saccades;
        end
    end
    end
end
end
end
% Neural_Properties
for counter_pCell = 1 : 1 : num_pCells
    ALL_PCELL_all_saccades.Neural_Properties.SS_num(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_num;
    ALL_PCELL_all_saccades.Neural_Properties.SS_duration(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_duration;
    ALL_PCELL_all_saccades.Neural_Properties.SS_firing_rate(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_firing_rate;
%     ALL_PCELL_all_saccades.Neural_Properties.SS_time(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_time;
    ALL_PCELL_all_saccades.Neural_Properties.SS_waveform(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.SS_waveform;
    ALL_PCELL_all_saccades.Neural_Properties.CS_num(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_num;
    ALL_PCELL_all_saccades.Neural_Properties.CS_firing_rate(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_firing_rate;
%     ALL_PCELL_all_saccades.Neural_Properties.CS_time(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_time;
    ALL_PCELL_all_saccades.Neural_Properties.CS_waveform(counter_pCell,:) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.CS_waveform;
    field_names_Corr_data = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.Corr_data);
    for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
        field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
        data_field_name_Corr_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties.Corr_data.(field_name_Corr_data);
        ALL_PCELL_all_saccades.Neural_Properties.Corr_data.(field_name_Corr_data)(counter_pCell,:) = data_field_name_Corr_data;
    end
end
fprintf(' --> Completed. \n')
end

%% function avg_over_amplitude
function ALL_PCELL_all_saccades_ang = avg_over_amplitude(ALL_PCELL_all_saccades)
%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades_ang', ' ... ']);
clearvars ALL_PCELL_all_saccades_ang
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 1);
span_width = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);
    
    for counter_ang = 1 : num_ang_bin
        train_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){1, counter_ang}));
        velocity_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){1, counter_ang}));
        num_saccades_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang}));
        if size(num_saccades_all_ang, 2) == 1
            num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
        end
        for counter_amp = 1 : num_amp_bin
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
            if size(num_saccades, 2) == 1
                num_saccades = repmat(num_saccades, 1, span_width);
            end
            train_data_all_ang = train_data_all_ang + (train_data .* num_saccades);
            velocity_data_all_ang = velocity_data_all_ang + (velocity_data .* num_saccades);
            num_saccades_all_ang = num_saccades_all_ang + num_saccades;
        end
        train_data_all_ang = train_data_all_ang ./ num_saccades_all_ang;
        train_data_all_ang(isnan(train_data_all_ang)) = 0;
        velocity_data_all_ang = velocity_data_all_ang ./ num_saccades_all_ang;
        velocity_data_all_ang(isnan(velocity_data_all_ang)) = 0;
        ALL_PCELL_all_saccades_ang.(variable_name).([spikeType_name '_train_' indType_name]){1, counter_ang} = train_data_all_ang;
        ALL_PCELL_all_saccades_ang.(variable_name).([spikeType_name '_velocity_' indType_name]){1, counter_ang} = velocity_data_all_ang;
        ALL_PCELL_all_saccades_ang.(variable_name).([spikeType_name '_num_sac_' indType_name]){1, counter_ang} = num_saccades_all_ang;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%% function avg_over_angle
function ALL_PCELL_all_saccades_amp = avg_over_angle(ALL_PCELL_all_saccades)
%% Loop over pCells
fprintf(['Building ALL_PCELL_all_saccades_amp', ' ... ']);
clearvars ALL_PCELL_all_saccades_amp
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
num_pCells = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 1);
span_width = size(ALL_PCELL_all_saccades.primSac.SS_train_start{1,1}, 2);
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);
    
    for counter_amp = 1 : num_amp_bin
        train_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, 1}));
        velocity_data_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, 1}));
        num_saccades_all_ang = zeros(size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1}));
        if size(num_saccades_all_ang, 2) == 1
            num_saccades_all_ang = repmat(num_saccades_all_ang, 1, span_width);
        end
        for counter_ang = 1 : num_ang_bin
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang};
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang};
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang};
            if size(num_saccades, 2) == 1
                num_saccades = repmat(num_saccades, 1, span_width);
            end
            train_data_all_ang = train_data_all_ang + (train_data .* num_saccades);
            velocity_data_all_ang = velocity_data_all_ang + (velocity_data .* num_saccades);
            num_saccades_all_ang = num_saccades_all_ang + num_saccades;
        end
        train_data_all_ang = train_data_all_ang ./ num_saccades_all_ang;
        train_data_all_ang(isnan(train_data_all_ang)) = 0;
        velocity_data_all_ang = velocity_data_all_ang ./ num_saccades_all_ang;
        velocity_data_all_ang(isnan(velocity_data_all_ang)) = 0;
        ALL_PCELL_all_saccades_amp.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, 1} = train_data_all_ang;
        ALL_PCELL_all_saccades_amp.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, 1} = velocity_data_all_ang;
        ALL_PCELL_all_saccades_amp.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, 1} = num_saccades_all_ang;
    end
end
end
end
fprintf(' --> Completed. \n')
end

%% function build_ALL_PCELL_all_saccades_tuned
function ALL_PCELL_all_saccades_tuned = build_ALL_PCELL_all_saccades_tuned(ALL_PCELL_all_saccades)
%% find CS-on for each cell
ALL_PCELL_all_saccades_ang = avg_over_amplitude(ALL_PCELL_all_saccades);
num_ang_bin = size(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start,2);
num_pCells = size(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start{1,1},1);
% METHOD 1
%{
CS_prob_tgtCuePres = zeros(num_pCells, num_ang_bin);
CS_prob_tgtPrimSac = zeros(num_pCells, num_ang_bin);
CS_prob_tgtEndPres = zeros(num_pCells, num_ang_bin);
CS_prob_tgtCorrSac = zeros(num_pCells, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_tgtCuePres(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtCuePres.CS_train_start{1,counter_ang}(:,300:500), 2);
    CS_prob_tgtPrimSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtPrimSac.CS_train_start{1,counter_ang}(:,100:300), 2);
    CS_prob_tgtEndPres(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtPrimSac.CS_train_finish{1,counter_ang}(:,300:500), 2);
    CS_prob_tgtCorrSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.tgtCorrSac.CS_train_start{1,counter_ang}(:,100:300), 2);
end
CS_prob_tgtCuePres(1:51, [2 4 6 8]) = nan;
CS_prob_tgtPrimSac(1:51, [2 4 6 8]) = nan;
CS_prob_tgtEndPres(1:51, [2 4 6 8]) = nan;
CS_prob_tgtCorrSac(1:51, [2 4 6 8]) = nan;
CS_prob_avg = (CS_prob_tgtCuePres + CS_prob_tgtPrimSac + CS_prob_tgtEndPres + CS_prob_tgtCorrSac) ./ 4;
hFig = figure(2); clf(hFig); hold('on');
plot(nanmean(CS_prob_tgtCuePres))
plot(nanmean(CS_prob_tgtPrimSac))
plot(nanmean(CS_prob_tgtEndPres))
plot(nanmean(CS_prob_tgtCorrSac))
plot(nanmean(CS_prob_avg), 'k')
ylim([0.14 0.24])
%}
% METHOD 2 
%
CS_prob_primSac = zeros(num_pCells, num_ang_bin);
CS_prob_corrSac = zeros(num_pCells, num_ang_bin);
for counter_ang = 1 : num_ang_bin
    CS_prob_primSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.primSac.CS_train_start{1,counter_ang}(:,100:300), 2);
    CS_prob_corrSac(:,counter_ang) = nansum(ALL_PCELL_all_saccades_ang.corrSac.CS_train_start{1,counter_ang}(:,100:300), 2);
end
CS_prob_primSac(1:51, [2 4 6 8]) = nan;
CS_prob_corrSac(1:51, [2 4 6 8]) = nan;
CS_prob_avg = (CS_prob_primSac + CS_prob_corrSac) ./ 2;
hFig = figure(1); clf(hFig); hold('on');
plot(nanmean(CS_prob_primSac))
plot(nanmean(CS_prob_corrSac))
plot(nanmean(CS_prob_avg), 'k')
ylim([0.14 0.24])
%}
[~, idx_CS_max] = nanmax(CS_prob_avg, [], 2);
idx_CS_tuning = zeros(num_pCells, num_ang_bin);
for counter_pCell = 1 : 1 : num_pCells
    CS_on_index = idx_CS_max(counter_pCell);
    CS_on_index = CS_on_index - 1; % CS_on_index should be in 0-index format
    if (CS_on_index == 8); CS_on_index = 0; end
    idx_CS_tuning(counter_pCell, :) = mod((CS_on_index : 1 : CS_on_index+7), 8) + 1;
end
%% build ALL_PCELL_all_saccades_tuned
fprintf(['Building ALL_PCELL_all_saccades_tuned', ' ... ']);
variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask', 'tgtCuePres', 'tgtPrimSac', 'tgtCorrSac'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};
clearvars ALL_PCELL_all_saccades_tuned
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    num_amp_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),1);
    num_ang_bin = size(ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]),2);

    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        for counter_pCell = 1 : 1 : num_pCells
            idx_ang_ = idx_CS_tuning(counter_pCell, counter_ang);
            
            train_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, idx_ang_}(counter_pCell,:);
            velocity_data = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, idx_ang_}(counter_pCell,:);
            num_saccades = ALL_PCELL_all_saccades.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, idx_ang_}(counter_pCell,:);
            
            ALL_PCELL_all_saccades_tuned.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = train_data;
            ALL_PCELL_all_saccades_tuned.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = velocity_data;
            ALL_PCELL_all_saccades_tuned.(variable_name).([spikeType_name '_num_sac_' indType_name]){counter_amp, counter_ang}(counter_pCell,:) = num_saccades;
        end
    end
    end
end
end
end
ALL_PCELL_all_saccades_tuned.idx_CS_tuning = idx_CS_tuning;
fprintf(' --> Completed. \n')

end

%% function build_TRAIN_DATA
function [TRAIN_DATA] = build_TRAIN_DATA(file_path, file_name)
%% Handle inputs
if nargin < 1
    [file_name_,file_path] = uigetfile([pwd filesep '*.psort'], 'Select psort file');
    [~,file_name,~] = fileparts(file_name_);
end

%% load EPHYS sorted DATA
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
file_name = [file_name '.psort'];
fprintf(['Loading ', file_name, ' ... ']);
% EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load EPHYS EVENT DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_aligned.mat'];
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.EPHYS_time_1K  = reshape(EPHYS.CH_EVE.EPHYS_time_1K ,[], 1);
EPHYS.CH_EVE.BEHAVE_time_1K = reshape(EPHYS.CH_EVE.BEHAVE_time_1K,[], 1);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_ANALYZED.mat'];
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path file_name]);
fprintf(' --> Completed. \n')

%% build EPHYS.CH_sorted from DATA_PSORT
ch_data = double(DATA_PSORT.topLevel_data.ch_data);
ch_time = double(DATA_PSORT.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT.topLevel_data.cs_index)));
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

waveform_inds_span = ((-60+1) : 1 : (120));
SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
CS_waveform = ch_data(CS_inds);

SS_waveform  = reshape(SS_waveform,[], length(waveform_inds_span));
CS_waveform  = reshape(CS_waveform,[], length(waveform_inds_span));

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data.CS_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_sorted.Corr_data.CS_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_sorted.Corr_data.SS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO    = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO    = Corr_data.CS_CSxSS_AUTO;
fprintf(' --> Completed. \n')

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K; % EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K; % 
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;

fprintf(' --> Completed. \n')

%% Building TRAIN_DATA
clearvars -except EPHYS BEHAVE TRAIN_DATA
fprintf(['Building TRAIN_DATA', ' ... ']);
TRAIN_DATA = struct;
inds_span    = ((-300+1) : 1 : (300))';
amp_edges = [-0.5 2 4 6 8 10 50];
ang_edges = -202.5 : 45 : +202.5;
TRAIN_DATA.inds_span = inds_span;
TRAIN_DATA.amp_edges = amp_edges;
TRAIN_DATA.ang_edges = ang_edges;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;

ang_edges_last_bin_id = length(ang_edges) - 1;

[~, ~, amp_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_amp_m, amp_edges);
[~, ~, ang_bin] = histcounts(BEHAVE.SACS_ALL_DATA.eye_r_ang, ang_edges);
ang_bin(ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

variable_list = {'primSac', 'corrSac', 'toTgt', 'fromTgt', 'lowAmp', 'toTgtStr', 'fromCenter', 'aroundCenter', 'nonTask'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};


for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = BEHAVE.SACS_ALL_DATA.(['is_' variable_name]);
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
    ind_indType_name = BEHAVE.SACS_ALL_DATA.(['ind_' indType_name]);
    ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (ang_bin == counter_ang);
        is_amp_bin   = (amp_bin == counter_amp);
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(1, length(inds_span));
            velocity_data = zeros(1, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = nanmean( spike_train_data(inds_train) );
        velocity_data = nanmean( velocity_trace_data(inds_velocity) );
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

%% Add tgtCuePres
clearvars -except EPHYS BEHAVE TRAIN_DATA
inds_span    = TRAIN_DATA.inds_span;
amp_edges = TRAIN_DATA.amp_edges;
ang_edges = TRAIN_DATA.ang_edges;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;
ang_edges_last_bin_id = length(ang_edges) - 1;

range_trials = BEHAVE.aligned.BEHAVE_range_trials;
tgt_amp_m = sqrt( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)).^2 + ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)).^2 ...
    );
tgt_ang = atan2d( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)) , ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)) ...
    );

[~, ~, amp_bin] = histcounts(tgt_amp_m, amp_edges);
[~, ~, ang_bin] = histcounts(tgt_ang  , ang_edges);
ang_bin(ang_bin == ang_edges_last_bin_id) = 1;

num_amp_bin = length(amp_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

variable_list = {'tgtCuePres'};
indType_list = {'start', 'vmax', 'finish'};
spikeType_list = {'SS', 'CS'};

% TRAIN_DATA = struct;
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    indType_name = indType_list{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = reshape(BEHAVE.aligned.BEHAVE_validity_prim, 1, []);
    ind_indType_name = BEHAVE.aligned.BEHAVE_ind_cue_present;
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
%     ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    ind_converted = EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (ang_bin == counter_ang);
        is_amp_bin   = (amp_bin == counter_amp);
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(1, length(inds_span));
            velocity_data = zeros(1, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = nanmean( spike_train_data(inds_train) );
        velocity_data = nanmean( velocity_trace_data(inds_velocity) );
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

%% Add tgtPrimSac tgtCorrSac
clearvars -except EPHYS BEHAVE TRAIN_DATA
inds_span    = TRAIN_DATA.inds_span;
amp_edges = TRAIN_DATA.amp_edges;
ang_edges = TRAIN_DATA.ang_edges;
length_train_data_ = length(EPHYS.CH_EVE.EPHYS_time_1K);
length_velocity_data_ = length(BEHAVE.aligned.time_1K);
velocity_trace_data = BEHAVE.aligned.eye_r_vm_filt;
ang_edges_last_bin_id = length(ang_edges) - 1;

range_trials = BEHAVE.aligned.BEHAVE_range_trials;
tgtPrimSac_amp_m = sqrt( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)).^2 + ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)).^2 ...
    );
tgtPrimSac_ang = atan2d( ...
    (BEHAVE.TRIALS_DATA.cue_y(range_trials) - BEHAVE.TRIALS_DATA.start_y(range_trials)) , ...
    (BEHAVE.TRIALS_DATA.cue_x(range_trials) - BEHAVE.TRIALS_DATA.start_x(range_trials)) ...
    );

tgtCorrSac_amp_m = sqrt( ...
    (BEHAVE.TRIALS_DATA.end_y(range_trials) - BEHAVE.TRIALS_DATA.cue_y(range_trials)).^2 + ...
    (BEHAVE.TRIALS_DATA.end_x(range_trials) - BEHAVE.TRIALS_DATA.cue_x(range_trials)).^2 ...
    );
tgtCorrSac_ang = atan2d( ...
    (BEHAVE.TRIALS_DATA.end_y(range_trials) - BEHAVE.TRIALS_DATA.cue_y(range_trials)) , ...
    (BEHAVE.TRIALS_DATA.end_x(range_trials) - BEHAVE.TRIALS_DATA.cue_x(range_trials)) ...
    );

[~, ~, tgtPrimSac_amp_m_bin] = histcounts(tgtPrimSac_amp_m, amp_edges);
[~, ~, tgtPrimSac_ang_bin] = histcounts(tgtPrimSac_ang  , ang_edges);
tgtPrimSac_ang_bin(tgtPrimSac_ang_bin == ang_edges_last_bin_id) = 1;

[~, ~, tgtCorrSac_amp_m_bin] = histcounts(tgtCorrSac_amp_m, amp_edges);
[~, ~, tgtCorrSac_ang_bin] = histcounts(tgtCorrSac_ang  , ang_edges);
tgtCorrSac_ang_bin(tgtCorrSac_ang_bin == ang_edges_last_bin_id) = 1;

bin_data.tgtPrimSac_amp_m_bin = tgtPrimSac_amp_m_bin;
bin_data.tgtPrimSac_ang_bin   = tgtPrimSac_ang_bin;
bin_data.tgtCorrSac_amp_m_bin = tgtCorrSac_amp_m_bin;
bin_data.tgtCorrSac_ang_bin   = tgtCorrSac_ang_bin;

num_amp_bin = length(amp_edges) - 1;
num_ang_bin = length(ang_edges) - 2;

variable_list = {'tgtPrimSac', 'tgtCorrSac'};
variable_list2 = {'primSac', 'corrSac'};
indType_list = {'start', 'vmax', 'finish'};
indType_list2 = {'onset', 'vmax', 'offset'};
spikeType_list = {'SS', 'CS'};

BEHAVE.aligned.BEHAVE_validity_primSac = BEHAVE.aligned.BEHAVE_validity_prim;
BEHAVE.aligned.BEHAVE_validity_corrSac = BEHAVE.aligned.BEHAVE_validity_corr;% & BEHAVE.aligned.BEHAVE_validity_prim;
% TRAIN_DATA = struct;
for counter_variable  = 1 : length(variable_list)
for counter_indType   = 1 : length(indType_list)
for counter_spikeType = 1 : length(spikeType_list)
    variable_name = variable_list{counter_variable};
    variable_name2 = variable_list2{counter_variable};
    indType_name = indType_list{counter_indType};
    indType_name2 = indType_list2{counter_indType};
    spikeType_name = spikeType_list{counter_spikeType};
    
    TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]) = cell(num_amp_bin, num_ang_bin);
    TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name]) = zeros(num_amp_bin, num_ang_bin);
    is_variable_name = reshape(BEHAVE.aligned.(['BEHAVE_validity_' variable_name2]), 1, []);
    ind_indType_name = BEHAVE.aligned.(['BEHAVE_ind_' variable_name2 '_' indType_name2]);
    spike_train_data = EPHYS.CH_EVE.(['EPHYS_' spikeType_name '_train_1K']);
    ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(ind_indType_name);
    
    for counter_amp = 1 : num_amp_bin
    for counter_ang = 1 : num_ang_bin
        is_ang_bin   = (bin_data.([variable_name '_ang_bin'])  == counter_ang );
        is_amp_bin   = (bin_data.([variable_name '_amp_m_bin'])  == counter_amp );
        if contains(variable_name,'tgtPrimSac') && contains(indType_name,'finish')
            % use tgtCorrSac angle for primSac_finish
            is_ang_bin   = (bin_data.(['tgtCorrSac' '_ang_bin'])  == counter_ang );
            is_amp_bin   = (bin_data.(['tgtPrimSac' '_amp_m_bin'])  == counter_amp );
        end
        ind_train    = ind_converted(is_variable_name & is_ang_bin & is_amp_bin);
        ind_velocity = ind_indType_name(is_variable_name & is_ang_bin & is_amp_bin);
        if isempty(ind_train)
            train_data = zeros(1, length(inds_span));
            velocity_data = zeros(1, length(inds_span));
            num_saccades = 0;
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        inds_train = repmat( ind_train(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_train), 1);
        inds_train( inds_train < 1 ) = 1;
        inds_train( inds_train > length_train_data_ ) = length_train_data_;
        inds_velocity = repmat( ind_velocity(:), 1, length(inds_span)) + repmat(inds_span(:)', length(ind_velocity), 1);
        inds_velocity( inds_velocity < 1 ) = 1;
        inds_velocity( inds_velocity > length_velocity_data_ ) = length_velocity_data_;
        if length(ind_train) == 1
            train_data = reshape( spike_train_data(inds_train) , 1, []);
            velocity_data = reshape( velocity_trace_data(inds_velocity) , 1, []);
            num_saccades = length(ind_train);
            TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
            TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
            continue;
        end
        train_data = nanmean( spike_train_data(inds_train) );
        velocity_data = nanmean( velocity_trace_data(inds_velocity) );
        num_saccades = length(ind_train);
        TRAIN_DATA.(variable_name).([spikeType_name '_train_' indType_name]){counter_amp, counter_ang} = train_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_velocity_' indType_name]){counter_amp, counter_ang} = velocity_data;
        TRAIN_DATA.(variable_name).([spikeType_name '_num_sac_' indType_name])(counter_amp, counter_ang) = num_saccades;
    end
    end
end
end
end

%% Add Neural_Properties
clearvars -except EPHYS BEHAVE TRAIN_DATA
TRAIN_DATA.Neural_Properties = struct();
if ~isempty(EPHYS.CH_sorted.SS_data.SS_time)
    TRAIN_DATA.Neural_Properties.SS_num = length(EPHYS.CH_sorted.SS_data.SS_time);
    TRAIN_DATA.Neural_Properties.SS_duration = (EPHYS.CH_sorted.SS_data.SS_time(end)) - (EPHYS.CH_sorted.SS_data.SS_time(1));
    TRAIN_DATA.Neural_Properties.SS_firing_rate = TRAIN_DATA.Neural_Properties.SS_num / TRAIN_DATA.Neural_Properties.SS_duration;
    TRAIN_DATA.Neural_Properties.SS_time = EPHYS.CH_sorted.SS_data.SS_time;
    TRAIN_DATA.Neural_Properties.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform;
else
    TRAIN_DATA.Neural_Properties.SS_num = 0;
    TRAIN_DATA.Neural_Properties.SS_duration = (EPHYS.CH_sorted.CS_data.CS_time(end)) - (EPHYS.CH_sorted.CS_data.CS_time(1));
    TRAIN_DATA.Neural_Properties.SS_firing_rate = 0;
    TRAIN_DATA.Neural_Properties.SS_time = [];
    TRAIN_DATA.Neural_Properties.SS_waveform = [];
end

if ~isempty(EPHYS.CH_sorted.CS_data.CS_time)
    TRAIN_DATA.Neural_Properties.CS_num = length(EPHYS.CH_sorted.CS_data.CS_time);
    TRAIN_DATA.Neural_Properties.CS_firing_rate = TRAIN_DATA.Neural_Properties.CS_num / TRAIN_DATA.Neural_Properties.SS_duration;
    TRAIN_DATA.Neural_Properties.CS_time = EPHYS.CH_sorted.CS_data.CS_time;
    TRAIN_DATA.Neural_Properties.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform;
else
    TRAIN_DATA.Neural_Properties.CS_num = 0;
    TRAIN_DATA.Neural_Properties.CS_firing_rate = 0;
    TRAIN_DATA.Neural_Properties.CS_time = [];
    TRAIN_DATA.Neural_Properties.CS_waveform = [];
end
TRAIN_DATA.Neural_Properties.Corr_data = EPHYS.CH_sorted.Corr_data;

%% Remove some fields
rmfields_list = {'inds_span', 'amp_edges', 'ang_edges'};
TRAIN_DATA = rmfield(TRAIN_DATA,rmfields_list);
fprintf(' --> Completed. \n')

end

%% function ESN_smooth
function smooth_data_ = ESN_smooth(data_, dim)
% smooth data using 2nd order Savitzky-Golay filter with 21 points
% if data_ is a matrix, the method will smooth each column by default or smooth along dim.
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
if nargin < 2
    dim = 1;
end
if size(data_, 1) == 1
    smooth_data_ = reshape(smooth(data_, 21, 'sgolay', 2), 1, []);
elseif size(data_, 2) == 1
    smooth_data_ = reshape(smooth(data_, 21, 'sgolay', 2), [], 1);
else
    smooth_data_ = nan(size(data_));
    if dim == 1
        % smooth columns
        for counter_col = 1 : size(data_, 2)
            smooth_data_(:, counter_col) = reshape(smooth(data_(:, counter_col), 21, 'sgolay', 2), [], 1);
        end
    elseif dim == 2
        % smooth rows
        for counter_row = 1 : size(data_, 1)
            smooth_data_(counter_row, :) = reshape(smooth(data_(counter_row, :), 21, 'sgolay', 2), 1, []);
        end
    end
    
end
end

%% function ESN_correlogram
function Corr_data = ESN_correlogram(SS_time, CS_time)
bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';

if (~isempty(CS_time)) && (~isempty(SS_time))
    ch_time_min = min([SS_time(1) CS_time(1)]);
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max([SS_time(end) CS_time(end)]) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(CS_time))
    ch_time_min = min(  CS_time(1) );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max(  CS_time(end) ) + 2.0;
    
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(SS_time))
    ch_time_min = min( SS_time(1)  );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max( SS_time(end)  ) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
end

% SSxSS_AUTO
if (~isempty(SS_time))
    CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
    % SSxSS correlogram
    SSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    ss_inds_span     = repmat(inds_span(:)',     size(SS_time(:),1), 1);
    ss_bin_size_time = repmat(bin_size_time(:)', size(SS_time(:),1), 1);
else
    SSxSS_AUTO       = false(0, length(inds_span(:)'));
    ss_inds_span     = nan(0, length(inds_span(:)'));
    ss_bin_size_time = nan(0, 1);
end

% CSxSS_WITHIN
if (~isempty(CS_time)) && (~isempty(SS_time))
    CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
    % CSxSS correlogram
    CSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    cs_inds_span     = repmat(inds_span(:)',     size(CS_time(:),1), 1);
    cs_bin_size_time = repmat(bin_size_time(:)', size(CS_time(:),1), 1);
else
    CSxSS_AUTO       = false(0, length(inds_span(:)'));
    cs_inds_span     = nan(0, length(inds_span(:)'));
    cs_bin_size_time = nan(0, 1);
end

Corr_data = struct;
Corr_data.CS_inds_span     = cs_inds_span;
Corr_data.CS_bin_size_time = cs_bin_size_time;
Corr_data.SS_inds_span     = ss_inds_span;
Corr_data.SS_bin_size_time = ss_bin_size_time;
Corr_data.SS_SSxSS_AUTO    = SSxSS_AUTO;
Corr_data.CS_CSxSS_AUTO    = CSxSS_AUTO;
end

%% function RE-RUN ESN_monkey_behavior_all_saccades
function re_run_ESN_monkey_behavior_all_saccades(pCell_list, path_data_monkey_sorted)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% RE-RUN ESN_monkey_behavior_all_saccades
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        file_name = [file_name_cell(1:13) '_ANALYZED.mat'];
        ESN_monkey_behavior_all_saccades(file_path, file_name);
    end
end

end

%% function scratch_plot1
function scratch_plot1
%%
hFig = figure(1);
clf(hFig)
subplot(1, 2, 1)
hold on
variable_name1 =  'tgtPrimSac';
variable_name2 =  'CS_train_start';
variable_name2V = 'CS_velocity_start';
variable_name3 =  'CS_num_sac_start';
ind_ang = 1;
range_ = 001:600;
num_pCells = 110;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned_ang;

data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){1, ind_ang}*1000;% - ALL_PCELL_DATA_1.(variable_name1).(variable_name2){1, 5}*1000;
data_ = data_(:, range_);
plot(ESN_smooth(nanmean(data_)+(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot(ESN_smooth(nanmean(data_)-(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot(ESN_smooth(nanmean(data_)), 'k', 'LineWidth', 2)
% set(gca, 'XTick', 0:50:400)
% ylim([40 90])
ylim([0 3])

subplot(1, 2, 2)
hold on
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){1, ind_ang};% - ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){1, 5};
data_ = data_(:, range_);
plot((nanmean(data_)+(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot((nanmean(data_)-(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot((nanmean(data_)), 'k', 'LineWidth', 2)
% set(gca, 'XTick', 0:50:400)
ylim([0 600])

ESN_Beautify_Plot
end

%% function scratch_plot2
function scratch_plot2
%% 
hFig = figure(2);
clf(hFig)
subplot(1, 2, 1)
hold on
variable_name1 =  'toTgtStr';
variable_name2 =  'SS_train_start';
variable_name2V = 'SS_velocity_start';
variable_name3 =  'SS_num_sac_start';
ind_ang = 5;
range_ = 201:400;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned;
ALL_PCELL_DATA_2 = ALL_PCELL_all_saccades_tuned_ang;

data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){2, ind_ang}(:,range_)*1000;
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){3, ind_ang}(:,range_)*1000;
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){4, ind_ang}(:,range_)*1000;
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_)))
% data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){5, ind_ang}(:,range_)*1000;
% data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang}==0), :) = nan;
% plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_2.(variable_name1).(variable_name2){1, ind_ang}(:,range_)*1000;
plot(ESN_smooth(nanmean(data_)), 'k', 'LineWidth', 2)
ylim([50 80])

subplot(1, 2, 2)
hold on
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){2, ind_ang}(:,range_);
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang}(:,1)==0), :) = nan;
plot((nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){3, ind_ang}(:,range_);
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang}(:,1)==0), :) = nan;
plot((nanmean(data_)))
data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){4, ind_ang}(:,range_);
data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang}(:,1)==0), :) = nan;
plot((nanmean(data_)))
% data_ = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){5, ind_ang}(:,range_);
% data_((ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang}==0), :) = nan;
% plot(ESN_smooth(nanmean(data_)))
data_ = ALL_PCELL_DATA_2.(variable_name1).(variable_name2V){1, ind_ang}(:,range_);
plot((nanmean(data_)), 'k', 'LineWidth', 2)
ylim([0 600])

ESN_Beautify_Plot
end

%% function scratch_plot3
function scratch_plot3
%% 
hFig = figure(3);
clf(hFig)
subplot(1, 2, 1)
hold on
variable_name1 =  'toTgtStr';
variable_name2 =  'SS_train_start';
variable_name2V = 'SS_velocity_start';
variable_name3 =  'SS_num_sac_start';
ind_ang_1 = 5;
ind_ang_2 = 1;
range_ = 201:400;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned;
ALL_PCELL_DATA_2 = ALL_PCELL_all_saccades_tuned_ang;

data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){2, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){2, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_2}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_1-data_2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){3, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){3, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_2}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_1-data_2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){4, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){4, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_2}(:,1)==0), :) = nan;
plot(ESN_smooth(nanmean(data_1-data_2)))
% data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){5, ind_ang_1}(:,range_);
% data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_1}(:,1)==0), :) = nan;
% data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2){5, ind_ang_2}(:,range_);
% data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_2}(:,1)==0), :) = nan;
% plot(ESN_smooth(nanmean(data_1-data_2)))
data_1 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2){1, ind_ang_1}(:,range_);
data_2 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2){1, ind_ang_2}(:,range_);
plot(ESN_smooth(nanmean(data_1-data_2)), 'k')

subplot(1, 2, 2)
hold on
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){2, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){2, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){2, ind_ang_2}(:,1)==0), :) = nan;
plot((nanmean((data_1+data_2)./2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){3, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){3, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){3, ind_ang_2}(:,1)==0), :) = nan;
plot((nanmean((data_1+data_2)./2)))
data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){4, ind_ang_1}(:,range_);
data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_1}(:,1)==0), :) = nan;
data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){4, ind_ang_2}(:,range_);
data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){4, ind_ang_2}(:,1)==0), :) = nan;
plot((nanmean((data_1+data_2)./2)))
% data_1 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){5, ind_ang_1}(:,range_);
% data_1(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_1}(:,1)==0), :) = nan;
% data_2 = ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){5, ind_ang_2}(:,range_);
% data_2(( ALL_PCELL_DATA_1.(variable_name1).(variable_name3){5, ind_ang_2}(:,1)==0), :) = nan;
% plot((nanmean((data_1+data_2)./2)))
data_1 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2V){1, ind_ang_1}(:,range_);
data_2 = ALL_PCELL_DATA_2.(variable_name1).(variable_name2V){1, ind_ang_2}(:,range_);
plot((nanmean((data_1+data_2)./2)), 'k')
end


%% function scratch_plot4
function scratch_plot4
%%
hFig = figure(4);
clf(hFig)

variable_name_category = 'nonTask';
event_type = 'vmax';
spike_type = 'SS';
variable_name_spike    = [spike_type '_train_' event_type];
variable_name_velocity = [spike_type '_velocity_' event_type];
variable_name_num_sac  = [spike_type '_num_sac_' event_type];
variable_name_firing   = [spike_type '_firing_rate'];
ind_ang = 1;
ind_amp = 2;
range_ = 201:400;
num_pCells = 110;

ALL_PCELL_DATA_1 = ALL_PCELL_all_saccades_tuned;

num_sac = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_num_sac){ind_amp, ind_ang};
firing_rate = ALL_PCELL_all_saccades.Neural_Properties.(variable_name_firing);
train_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_spike){ind_amp, ind_ang}*1000;
train_data = train_data - repmat(firing_rate(:,1), 1, size(train_data, 2));
num_sac_2 = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_num_sac){ind_amp, 5};
train_data_2 = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_spike){ind_amp, 5}*1000;
train_data_2 = train_data_2 - repmat(firing_rate(:,1), 1, size(train_data_2, 2));

num_perm = 5000;
num_samples = num_pCells ; % 50; % 
length_span = size(train_data, 2);
train_data_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_ = train_data(inds_perm, :);
    train_data_(isnan(train_data_)) = 0;
    num_sac_ = num_sac(inds_perm, 1);
    firing_mean_ = num_sac_' * train_data_ ./ nansum(num_sac_);
    train_data_perm(counter_perm, :) = firing_mean_;
    
    inds_perm_2 = inds_perm; % randi(num_pCells,[num_samples 1]); % 
    train_data_2_ = train_data_2(inds_perm_2, :);
    train_data_2_(isnan(train_data_2_)) = 0;
    num_sac_2_ = num_sac_2(inds_perm_2, 1);
    firing_mean_2_ = num_sac_2_' * train_data_2_ ./ nansum(num_sac_2_);
    train_data_perm(counter_perm, :) = firing_mean_ - firing_mean_2_;
    
end
data_mean = ESN_smooth( nanmean(train_data_perm, 1) );
data_stdv = ESN_smooth( nanstd( train_data_perm, 0, 1) ) ./ sqrt(num_samples) * 3; % 95% CI
data_mean = data_mean(:, range_);
data_stdv = data_stdv(:, range_);

subplot(1, 2, 1)
hold on
y_axis_data_stdv = [(data_mean + data_stdv) flip(data_mean - data_stdv)];
x_axis_data_stdv = [(1:1:length(data_mean)) (length(data_mean):-1:1)];
plot(x_axis_data_stdv, (y_axis_data_stdv), 'k', 'LineWidth', 0.25)
% plot((data_mean - data_stdv), 'k', 'LineWidth', 0.25)
plot((data_mean), 'k', 'LineWidth', 1)
% set(gca, 'XTick', 0:200:600)
ylim([-20 30])
% ylim([-1 1.5])

subplot(1, 2, 2)
hold on
velocity_data = ALL_PCELL_DATA_1.(variable_name_category).(variable_name_velocity){ind_amp, ind_ang};% - ALL_PCELL_DATA_1.(variable_name1).(variable_name2V){1, 5};
velocity_data((num_sac(:,1)==0), :) = nan;
velocity_data = velocity_data(:, range_);
% plot((nanmean(data_)+(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
% plot((nanmean(data_)-(nanstd(data_)./sqrt(num_pCells))), 'k', 'LineWidth', 1)
plot((nanmean(velocity_data)), 'k', 'LineWidth', 1)
% set(gca, 'XTick', 0:200:600)
ylim([0 600])

ESN_Beautify_Plot(hFig, [3 1.5])

end

%%
function classify_pCells
%%
data_ = ALL_PCELL_all_saccades_tuned_ang_amp.primSac.SS_train_start{1,1};
SS_baseline = mean(data_(:,151:250), 2);
data_norm = data_ ./ repmat(SS_baseline, 1, size(data_, 2));
data_smooth = ESN_smooth(data_norm, 2);

[~, pca_mat, ~] = pca(data_smooth(:,250:400));
[reduction, umap, clusterIdentifiers, extras]=run_umap(data_smooth(:,250:400));

hFig = figure(1);
clf(hFig);
subplot(1, 3, 1)
plot(data_smooth')
ylim([0 inf])
subplot(1, 3, 2)
plot(pca_mat(:, 1), pca_mat(:, 2), '.k')

subplot(1, 3, 3)
plot(reduction(:, 1), reduction(:, 2), '.k')
%%
label_1 = reduction(:, 1) < 9.5;
hFig = figure(2);
clf(hFig);
subplot(1, 2, 1)
hold on
plot((-49:100)',data_smooth(label_1,251:400)', 'b')
plot((-49:100)',data_smooth(~label_1,251:400)', 'r')
xlabel('Saccade onset (ms)')
ylabel('Normalized SS firing rate')
ylim([0 inf])

subplot(1, 2, 2)
hold on
plot(reduction(label_1, 1), reduction(label_1, 2), 'ob')
plot(reduction(~label_1, 1), reduction(~label_1, 2), 'or')
xlabel('umap 1')
ylabel('umap 2')

ESN_Beautify_Plot

end