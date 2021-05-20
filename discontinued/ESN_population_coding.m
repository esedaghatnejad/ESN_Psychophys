%% function ESN_population_coding
function ESN_population_coding
%% build and save ALL_PCELL_COMPRESSED_DATA
% path to data_125d_sorted folder
% path_data_monkey_sorted = '/Volumes/ShadmehrLab/data_125d_sorted';
% path_data_monkey_sorted = '/Volumes/ShadmehrLab/data_59d_sorted';
% path_data_monkey_sorted = '/Volumes/ShadmehrLab/data_184d_sorted';
path_data_monkey_sorted = uigetdir;
% build pCell_list, this is a hard coded cell with the id of all of the pCells and the bundles

% pCell_list = build_pCell_list_Mirza_pre201906(); 
% pCell_list = build_pCell_list_Mirza_post201906();
% pCell_list = build_pCell_list_Ramon();

% pCell_list_1 = build_pCell_list_Mirza_pre201906(); pCell_list_2 = build_pCell_list_Mirza_post201906();
% pCell_list = vertcat(pCell_list_1, pCell_list_2);

% pCell_list_1 = build_pCell_list_Mirza_post201906(); pCell_list_2 = build_pCell_list_Ramon();
% pCell_list = vertcat(pCell_list_1, pCell_list_2);

pCell_list_1 = build_pCell_list_Mirza_pre201906(); 
pCell_list_2 = build_pCell_list_Mirza_post201906(); 
pCell_list_3 = build_pCell_list_Ramon();
pCell_list_4 = build_pCell_list_Mirza_post202011();
pCell_list = vertcat(pCell_list_1, pCell_list_2, pCell_list_3, pCell_list_4);
%%
% build ALL_PCELL_COMPRESSED_DATA, open plot data and put them together
ALL_PCELL_COMPRESSED_DATA = build_ALL_PCELL_COMPRESSED_DATA(pCell_list, path_data_monkey_sorted);

num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
path_data_monkey_sorted = [path_data_monkey_sorted filesep ALL_PCELL_name];
if ~exist(path_data_monkey_sorted, 'dir')
    mkdir(path_data_monkey_sorted);
end

% save ALL_PCELL_COMPRESSED_DATA
save([path_data_monkey_sorted, filesep, 'ALL_PCELL_COMPRESSED_DATA.mat'],...
    'ALL_PCELL_COMPRESSED_DATA',...
    '-v7.3')

%% load and plot ALL_PCELL_COMPRESSED_DATA
if ~exist('ALL_PCELL_COMPRESSED_DATA','var')
    % path to data_125d_sorted folder
    path_data_monkey_sorted = uigetdir;
    load([path_data_monkey_sorted, filesep, 'ALL_PCELL_COMPRESSED_DATA.mat'])
end
RASTER_DATA_ALL_PCELL_TUNED = build_RASTER_DATA_ALL_PCELL_TUNED(ALL_PCELL_COMPRESSED_DATA);
plot_RASTER_DATA_ALL_PCELL_TUNED(path_data_monkey_sorted, ALL_PCELL_COMPRESSED_DATA, RASTER_DATA_ALL_PCELL_TUNED);

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
pCell_list_Mirza_post201906( 1, 1:1) = {'190722_154305_03_sorted_RSh'};
pCell_list_Mirza_post201906( 2, 1:1) = {'190722_154305_07_sorted_PGH'};
pCell_list_Mirza_post201906( 3, 1:5) = {'190724_141024_04_sorted_PGH', '190724_142939_05_sorted_PGH', '190724_145556_06_sorted_PGH', '190724_151910_06_sorted_PGH', '190724_154019_06_sorted_PGH'};
pCell_list_Mirza_post201906( 4, 1:4) = {'190730_143843_01_sorted_PGH', '190730_150548_01_sorted_PGH', '190730_154007_01_sorted_PGH', '190730_161519_01_sorted_PGH'};
pCell_list_Mirza_post201906( 5, 1:1) = {'190801_144133_05_sorted_RSh'};
pCell_list_Mirza_post201906( 6, 1:1) = {'190801_151921_04_sorted_RSh'};
pCell_list_Mirza_post201906( 7, 1:1) = {'190801_161911_05_sorted_RSh'};
pCell_list_Mirza_post201906( 8, 1:2) = {'190812_140139_03_sorted_RSh', '190812_142255_03_sorted_RSh'};
pCell_list_Mirza_post201906( 9, 1:3) = {'190812_153354_02_sorted_RSh', '190812_160450_02_sorted_RSh', '190812_163649_02_sorted_RSh'};
pCell_list_Mirza_post201906(10, 1:3) = {'190815_161258_08_sorted_PGH', '190815_165738_08_sorted_PGH', '190815_171400_08_sorted_PGH'};

pCell_list_Mirza_post201906(11, 1:2) = {'190918_115126_04_sorted_RSh', '190918_122153_04_sorted_RSh'};
pCell_list_Mirza_post201906(12, 1:2) = {'190919_095155_02_sorted_RSh', '190919_100229_02_sorted_RSh'};
pCell_list_Mirza_post201906(13, 1:3) = {'190923_143725_02_sorted_PGH', '190923_150658_02_sorted_PGH', '190923_154446_02_sorted_PGH'};
pCell_list_Mirza_post201906(14, 1:3) = {'190925_145115_02_sorted_PGH', '190925_152228_02_sorted_PGH', '190925_155449_02_sorted_PGH'};
pCell_list_Mirza_post201906(15, 1:1) = {'190927_135905_01_sorted_RSh'};
pCell_list_Mirza_post201906(16, 1:1) = {'190927_135905_03_sorted_RSh'};

end

%% function build_pCell_list_Mirza_post202011
function pCell_list_Mirza_post202011 = build_pCell_list_Mirza_post202011()
%% build pCell_list
pCell_list_Mirza_post202011 = cell(0,10);

% pCell_list(0) = {'', '', '', '', '', '', '', '', '', ''};
pCell_list_Mirza_post202011( 1, 1:2) = {'201223_131532_07_sorted_RSh', '201223_134248_07_sorted_RSh'};
pCell_list_Mirza_post202011( 2, 1:1) = {'201223_131532_03_sorted_RSh'};
pCell_list_Mirza_post202011( 3, 1:2) = {'201223_134248_05_sorted_RSh', '201223_142804_05_sorted_RSh'};
pCell_list_Mirza_post202011( 4, 1:3) = {'201229_160007_04_sorted_RSh', '201229_161147_04_sorted_RSh', '201229_163338_04_sorted_RSh'};
pCell_list_Mirza_post202011( 5, 1:1) = {'201229_170800_03_sorted_RSh'};
pCell_list_Mirza_post202011( 6, 1:3) = {'201230_145929_04_sorted_RSh', '201230_151857_04_sorted_RSh', '201230_154916_04_sorted_RSh'};
pCell_list_Mirza_post202011( 7, 1:1) = {'201230_171821_04_sorted_RSh'};
pCell_list_Mirza_post202011( 8, 1:3) = {'210104_133940_04_sorted_RSh', '210104_135712_03_sorted_RSh', '210104_142745_03_sorted_RSh'};
pCell_list_Mirza_post202011( 9, 1:2) = {'210104_151914_01_sorted_RSh', '210104_153305_01_sorted_RSh'};
pCell_list_Mirza_post202011(10, 1:2) = {'210104_151914_02_sorted_RSh', '210104_153305_02_sorted_RSh'};

pCell_list_Mirza_post202011(11, 1:5) = {'210106_141642_03_sorted_RSh', '210106_142808_03_sorted_RSh', '210106_145826_03_sorted_RSh', '210106_150439_03_sorted_RSh', '210106_153603_03_sorted_RSh'};
pCell_list_Mirza_post202011(12, 1:5) = {'210106_141642_05_sorted_RSh', '210106_142808_05_sorted_RSh', '210106_145826_05_sorted_RSh', '210106_150439_05_sorted_RSh', '210106_153603_05_sorted_RSh'};
pCell_list_Mirza_post202011(13, 1:3) = {'210107_141037_04_sorted_RSh', '210107_143610_04_sorted_RSh', '210107_145924_04_sorted_RSh'};
pCell_list_Mirza_post202011(14, 1:1) = {'210107_154500_03_sorted_RSh'};
pCell_list_Mirza_post202011(15, 1:1) = {'210107_154500_06_sorted_RSh'};
pCell_list_Mirza_post202011(16, 1:1) = {'210107_154500_07_sorted_RSh'};
pCell_list_Mirza_post202011(17, 1:2) = {'210108_151142_06_sorted_RSh', '210108_153618_06_sorted_RSh'};
pCell_list_Mirza_post202011(18, 1:2) = {'210108_151142_07_sorted_RSh', '210108_153618_07_sorted_RSh'};
pCell_list_Mirza_post202011(19, 1:1) = {'210122_133338_04_sorted_RSh'};
pCell_list_Mirza_post202011(20, 1:5) = {'210122_142609_05_sorted_RSh', '210122_144322_05_sorted_RSh', '210122_145740_05_sorted_RSh', '210122_154050_05_sorted_RSh', '210122_161215_05_sorted_RSh'};

pCell_list_Mirza_post202011(21, 1:2) = {'210122_144322_06_sorted_RSh', '210122_145740_06_sorted_RSh'};
pCell_list_Mirza_post202011(22, 1:5) = {'210129_135158_04_sorted_RSh', '210129_142058_04_sorted_RSh', '210129_145659_04_sorted_RSh', '210129_151717_04_sorted_RSh', '210129_153949_04_sorted_RSh'};

end

%% function build_pCell_list_Ramon
function pCell_list_Ramon = build_pCell_list_Ramon()
%% build pCell_list
pCell_list_Ramon = cell(0,10);

% pCell_list(0) = {'', '', '', '', '', '', '', '', '', ''};
pCell_list_Ramon( 1, 1:3) = {'190829_131625_04_sorted_ESN', '190829_132447_04_sorted_ESN', '190829_133438_04_sorted_ESN'};
pCell_list_Ramon( 2, 1:2) = {'190830_114101_02_sorted_ESN', '190830_114850_02_sorted_ESN'};
pCell_list_Ramon( 3, 1:2) = {'190830_114101_04_sorted_ESN', '190830_114850_04_sorted_ESN'};
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
pCell_list_Ramon(32, 1:1) = {'191209_141244_03_sorted_JSP'};
pCell_list_Ramon(33, 1:1) = {'191209_141244_04_sorted_RSh'};
pCell_list_Ramon(34, 1:2) = {'200109_110904_04_sorted_RSh', '200109_112601_04_sorted_RSh'};
pCell_list_Ramon(35, 1:1) = {'200109_125642_03_sorted_RSh'};
pCell_list_Ramon(36, 1:3) = {'200122_125828_07_sorted_RSh', '200122_133020_07_sorted_RSh', '200122_135404_07_sorted_RSh'};
pCell_list_Ramon(37, 1:2) = {'200122_133020_04_sorted_RSh', '200122_135404_04_sorted_RSh'};
pCell_list_Ramon(38, 1:4) = {'200123_124808_07_sorted_PGH', '200123_131422_07_sorted_PGH', '200123_133835_07_sorted_PGH', '200123_140924_07_sorted_PGH'};
pCell_list_Ramon(39, 1:4) = {'200128_151251_04_sorted_RSh', '200128_151937_06_sorted_RSh', '200128_153045_06_sorted_RSh', '200128_154444_01_sorted_RSh'};
pCell_list_Ramon(40, 1:1) = {'200128_163856_04_sorted_RSh'};

pCell_list_Ramon(41, 1:4) = {'200129_094755_04_sorted_PGH', '200129_124934_04_sorted_PGH', '200129_131018_04_sorted_PGH', '200129_133230_04_sorted_PGH'};
pCell_list_Ramon(42, 1:3) = {'200204_162413_04_sorted_RSh', '200204_163427_04_sorted_RSh', '200204_164913_03_sorted_RSh'};
pCell_list_Ramon(43, 1:1) = {'200211_162022_03_sorted_RSh'};
pCell_list_Ramon(44, 1:5) = {'200303_153615_03_sorted_RSh', '200303_155436_03_sorted_RSh', '200303_163555_03_sorted_RSh', '200303_165850_03_sorted_RSh', '200303_172913_03_sorted_RSh'};
pCell_list_Ramon(45, 1:3) = {'200814_113452_05_sorted_RSh', '200814_114000_05_sorted_RSh', '200814_114425_05_sorted_RSh'};
pCell_list_Ramon(46, 1:3) = {'200914_122203_04_sorted_RSh', '200914_124954_04_sorted_RSh', '200914_132356_04_sorted_RSh'};
pCell_list_Ramon(47, 1:3) = {'200916_130045_02_sorted_PGH', '200916_133908_02_sorted_PGH', '200916_142429_02_sorted_PGH'};
pCell_list_Ramon(48, 1:1) = {'200923_133039_03_sorted_PGH'};
pCell_list_Ramon(49, 1:4) = {'201007_132359_04_sorted_PGH', '201007_133319_02_sorted_PGH', '201007_134353_02_sorted_PGH', '201007_135402_04_sorted_PGH'};

end

%% function build_ALL_PCELL_COMPRESSED_DATA
function ALL_PCELL_COMPRESSED_DATA = build_ALL_PCELL_COMPRESSED_DATA(pCell_list, path_data_monkey_sorted)
%% init vars
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
clearvars ALL_PCELL_COMPRESSED_DATA
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    plot_data_cell = struct();
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name = [pCell_list{counter_pCell, counter_recording} '_plot_data' '.mat']; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name(1:2);
        month_ = file_name(3:4);
        day_ = file_name(5:6);
        hour_ = file_name(8:9);
        minute_ = file_name(10:11);
        second_ = file_name(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_figs' filesep];
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        %% load plot_data
        plot_data_raw = load([file_path filesep file_name]);
        plot_data_raw.file_name = file_name;
        plot_data_raw.file_path = file_path;
        if length(plot_data_raw.file_name) == 41
            plot_data_raw.id          = plot_data_raw.file_name(1:16);
        elseif length(plot_data_raw.file_name) == 43
            plot_data_raw.id          = plot_data_raw.file_name(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        % accounting for a bug where 1 spike has vertical matrix instead of
        % horizontal
        if length(plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_ind) == 1
            plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_waveform = ...
                plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_waveform(:)';
        end
        if length(plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_ind) == 1
            plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_waveform = ...
                plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_waveform(:)';
        end
        %% concatenate dataset
        % concatenate file_name
        if ~isfield(plot_data_cell, 'file_name')
            plot_data_cell.id = cell(0,0);
            plot_data_cell.file_name = cell(0,0);
            plot_data_cell.file_path = cell(0,0);
        end
        plot_data_cell.id{counter_recording, 1}   = plot_data_raw.id;
        plot_data_cell.file_name{counter_recording, 1}   = plot_data_raw.file_name;
        plot_data_cell.file_path{counter_recording, 1}   = plot_data_raw.file_path;
        % concatenate raster_data
        field_names_plot_data_raw = fieldnames(plot_data_raw);
        field_names_plot_data_raw = field_names_plot_data_raw(contains(field_names_plot_data_raw, 'raster_data_'));
        for counter_field_names_plot_data_raw = 1 : 1 : length(field_names_plot_data_raw)
            field_name_plot_data_raw = field_names_plot_data_raw{counter_field_names_plot_data_raw};
            field_names_raster_data = fieldnames(plot_data_raw.(field_name_plot_data_raw));
            % the field does not exist in plot_data_cell
            if ~isfield(plot_data_cell, field_name_plot_data_raw)
                plot_data_cell.(field_name_plot_data_raw) = struct();
            end
            for counter_field_names_raster_data = 1 : 1 : length(field_names_raster_data)
                field_name_raster_data = field_names_raster_data{counter_field_names_raster_data};
                data_field_name_raster_data = plot_data_raw.(field_name_plot_data_raw).(field_name_raster_data);
                % the field does not exist in plot_data_cell
                if ~isfield(plot_data_cell.(field_name_plot_data_raw), field_name_raster_data)
                    plot_data_cell.(field_name_plot_data_raw).(field_name_raster_data) = [];
                end
                data_field_name_cell = plot_data_cell.(field_name_plot_data_raw).(field_name_raster_data);
                data_field_name_cell = vertcat(data_field_name_cell, data_field_name_raster_data);
                plot_data_cell.(field_name_plot_data_raw).(field_name_raster_data) = data_field_name_cell;
            end
            % add inds_span from plot_data to raster_data
            field_name_plot_data_event = ['plot' field_name_plot_data_raw(7:end)];
            field_name_inds_span = 'inds_span';
            data_field_name_raster_data = plot_data_raw.(field_name_plot_data_event).(field_name_inds_span);
            if ~isfield(plot_data_cell.(field_name_plot_data_raw), field_name_inds_span)
                plot_data_cell.(field_name_plot_data_raw).(field_name_inds_span) = [];
            end
            data_field_name_cell = plot_data_cell.(field_name_plot_data_raw).(field_name_inds_span);
            data_field_name_cell = vertcat(data_field_name_cell, data_field_name_raster_data);
            plot_data_cell.(field_name_plot_data_raw).(field_name_inds_span) = data_field_name_cell;
        end
        % concatenate Neural_Properties_data
        if ~isfield(plot_data_cell, 'Neural_Properties_data')
            plot_data_cell.Neural_Properties_data = struct();
        end
        if ~isfield(plot_data_cell.Neural_Properties_data, 'SS_num')
            plot_data_cell.Neural_Properties_data.SS_num = [];
            plot_data_cell.Neural_Properties_data.SS_duration = [];
            plot_data_cell.Neural_Properties_data.SS_firing_rate = [];
            plot_data_cell.Neural_Properties_data.SS_time = [];
            plot_data_cell.Neural_Properties_data.SS_waveform = [];
        end
        if ~isempty(plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time)
        plot_data_cell.Neural_Properties_data.SS_num(counter_recording, 1) = ...
            length(plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time);
        plot_data_cell.Neural_Properties_data.SS_duration(counter_recording, 1) = ...
            (plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time(end)) - (plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time(1));
        plot_data_cell.Neural_Properties_data.SS_firing_rate(counter_recording, 1) = ...
            plot_data_cell.Neural_Properties_data.SS_num(counter_recording, 1) / plot_data_cell.Neural_Properties_data.SS_duration(counter_recording, 1);
        SS_time_plot_data_raw = plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_time;
        SS_time_plot_data_cell = plot_data_cell.Neural_Properties_data.SS_time;
        SS_time_plot_data_cell = vertcat(SS_time_plot_data_cell, SS_time_plot_data_raw);
        plot_data_cell.Neural_Properties_data.SS_time = SS_time_plot_data_cell;
        SS_waveform_plot_data_raw = plot_data_raw.Neural_Properties_data.CH_sorted.SS_data.SS_waveform;
        SS_waveform_plot_data_cell = plot_data_cell.Neural_Properties_data.SS_waveform;
        SS_waveform_plot_data_cell = vertcat(SS_waveform_plot_data_cell, SS_waveform_plot_data_raw);
        plot_data_cell.Neural_Properties_data.SS_waveform = SS_waveform_plot_data_cell;
        else
        plot_data_cell.Neural_Properties_data.SS_num(counter_recording, 1) = 0;
        plot_data_cell.Neural_Properties_data.SS_duration(counter_recording, 1) = ...
            (plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time(end)) - (plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time(1));
        plot_data_cell.Neural_Properties_data.SS_firing_rate(counter_recording, 1) = 0;
        end
        
        if ~isfield(plot_data_cell.Neural_Properties_data, 'CS_num')
            plot_data_cell.Neural_Properties_data.CS_num = [];
            plot_data_cell.Neural_Properties_data.CS_firing_rate = [];
            plot_data_cell.Neural_Properties_data.CS_time = [];
            plot_data_cell.Neural_Properties_data.CS_waveform = [];
        end
        if ~isempty(plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time)
        plot_data_cell.Neural_Properties_data.CS_num(counter_recording, 1) = ...
            length(plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time);
        plot_data_cell.Neural_Properties_data.CS_firing_rate(counter_recording, 1) = ...
            plot_data_cell.Neural_Properties_data.CS_num(counter_recording, 1) / plot_data_cell.Neural_Properties_data.SS_duration(counter_recording, 1);
        CS_time_plot_data_raw = plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_time;
        CS_time_plot_data_cell = plot_data_cell.Neural_Properties_data.CS_time;
        CS_time_plot_data_cell = vertcat(CS_time_plot_data_cell, CS_time_plot_data_raw);
        plot_data_cell.Neural_Properties_data.CS_time = CS_time_plot_data_cell;
        CS_waveform_plot_data_raw = plot_data_raw.Neural_Properties_data.CH_sorted.CS_data.CS_waveform;
        CS_waveform_plot_data_cell = plot_data_cell.Neural_Properties_data.CS_waveform;
        CS_waveform_plot_data_cell = vertcat(CS_waveform_plot_data_cell, CS_waveform_plot_data_raw);
        plot_data_cell.Neural_Properties_data.CS_waveform = CS_waveform_plot_data_cell;
        else
        plot_data_cell.Neural_Properties_data.CS_num(counter_recording, 1) = 0;
        plot_data_cell.Neural_Properties_data.CS_firing_rate(counter_recording, 1) = 0;
        end
        % concatenate Neural_Properties_data.Corr_data
        if ~isfield(plot_data_cell.Neural_Properties_data, 'Corr_data')
            plot_data_cell.Neural_Properties_data.Corr_data = struct();
        end
        field_names_Corr_data = fieldnames(plot_data_raw.Neural_Properties_data.CH_sorted.Corr_data);
        for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
            field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
            data_field_name_Corr_data = plot_data_raw.Neural_Properties_data.CH_sorted.Corr_data.(field_name_Corr_data);
            if ~isfield(plot_data_cell.Neural_Properties_data.Corr_data, field_name_Corr_data)
                plot_data_cell.Neural_Properties_data.Corr_data.(field_name_Corr_data) = [];
            end
            data_field_name_cell = plot_data_cell.Neural_Properties_data.Corr_data.(field_name_Corr_data);
            data_field_name_cell = vertcat(data_field_name_cell, data_field_name_Corr_data);
            plot_data_cell.Neural_Properties_data.Corr_data.(field_name_Corr_data) = data_field_name_cell;
        end
    end
    %% Build plot_data_compress
    plot_data_compress = struct();
    plot_data_compress.id   = plot_data_cell.id;
    plot_data_compress.file_name   = plot_data_cell.file_name;
    plot_data_compress.file_path   = plot_data_cell.file_path;
    field_names_plot_data_cell = fieldnames(plot_data_cell);
    field_names_plot_data_cell = field_names_plot_data_cell(contains(field_names_plot_data_cell, 'raster_data_'));
    for counter_field_names_plot_data_cell = 1 : 1 : length(field_names_plot_data_cell)
        field_name_plot_data_cell = field_names_plot_data_cell{counter_field_names_plot_data_cell};
        field_names_raster_data = fieldnames(plot_data_cell.(field_name_plot_data_cell));
        for counter_field_names_raster_data = 1 : 1 : length(field_names_raster_data)
            field_name_raster_data = field_names_raster_data{counter_field_names_raster_data};
            data_field_name_raster_data = plot_data_cell.(field_name_plot_data_cell).(field_name_raster_data);
            plot_data_compress.(field_name_plot_data_cell).(field_name_raster_data) = data_field_name_raster_data;
            plot_data_compress.(field_name_plot_data_cell).(field_name_raster_data) = nanmean(data_field_name_raster_data, 1);
            if contains(field_name_raster_data, '_CS_')
                plot_data_compress.(field_name_plot_data_cell).([field_name_raster_data '_sparse']) = sparse(data_field_name_raster_data);
                plot_data_compress.(field_name_plot_data_cell).([field_name_raster_data '_numTrial']) = size(data_field_name_raster_data, 1);
            end
        end
    end
    
    plot_data_compress.Neural_Properties_data.SS_num = ...
        nansum(plot_data_cell.Neural_Properties_data.SS_num);
    plot_data_compress.Neural_Properties_data.SS_duration = ...
        nansum(plot_data_cell.Neural_Properties_data.SS_duration);
    plot_data_compress.Neural_Properties_data.SS_firing_rate = ...
        plot_data_compress.Neural_Properties_data.SS_num / plot_data_compress.Neural_Properties_data.SS_duration;
    plot_data_compress.Neural_Properties_data.SS_time = ...
        plot_data_cell.Neural_Properties_data.SS_time;
    plot_data_compress.Neural_Properties_data.SS_waveform = ...
        nanmean(plot_data_cell.Neural_Properties_data.SS_waveform, 1);
    
    plot_data_compress.Neural_Properties_data.CS_num = ...
        nansum(plot_data_cell.Neural_Properties_data.CS_num);
    plot_data_compress.Neural_Properties_data.CS_firing_rate = ...
        plot_data_compress.Neural_Properties_data.CS_num / plot_data_compress.Neural_Properties_data.SS_duration;
    plot_data_compress.Neural_Properties_data.CS_time = ...
        plot_data_cell.Neural_Properties_data.CS_time;
    plot_data_compress.Neural_Properties_data.CS_waveform = ...
        nanmean(plot_data_cell.Neural_Properties_data.CS_waveform, 1);
    
    field_names_Corr_data = fieldnames(plot_data_cell.Neural_Properties_data.Corr_data);
    for counter_field_names_Corr_data = 1 : 1 : length(field_names_Corr_data)
        field_name_Corr_data = field_names_Corr_data{counter_field_names_Corr_data};
        data_field_name_Corr_data = plot_data_cell.Neural_Properties_data.Corr_data.(field_name_Corr_data);
        plot_data_compress.Neural_Properties_data.Corr_data.(field_name_Corr_data) = nanmean(data_field_name_Corr_data, 1);
    end
    %% build plot_data_compress.CS_Tuning
    plot_data_compress.CS_Tuning = struct();
    event_names = {'cue_present', 'primSac_onset', 'primSac_offset', 'corrSac_onset'};
    range_inds_probability_list = [(301:500); (101:300); (301:500); (101:300);];
    direction_names = {'000', '045', '090', '135', '180', '225', '270', '315'};
    overall_prob = [];
    overall_numTrials = [];
    overall_vector_sum = [0, 0];
    for counter_event = 1 : 1 : length(event_names)
        event_name = event_names{counter_event};
        plot_data_compress.CS_Tuning.([event_name '_prob']) = [];
        event_vector_sum = [0, 0];
        event_numTrials = 0;
        for counter_direction = 1 : 1 : length(direction_names)
            direction_name = direction_names{counter_direction};
            range_inds_probability = range_inds_probability_list(counter_event, :);
            raster_data_ = full(plot_data_compress.(['raster_data_' event_name]).(['train_data_logic_CS_' direction_name '_sparse']));
            prob_ = nansum( nansum(raster_data_(:,range_inds_probability) ,2) > 0 ) / size(raster_data_, 1);
            plot_data_compress.CS_Tuning.([event_name '_' direction_name '_prob']) = prob_;
            data_event_prob = plot_data_compress.CS_Tuning.([event_name '_prob']);
            data_event_prob = horzcat(data_event_prob, prob_);
            plot_data_compress.CS_Tuning.([event_name '_prob']) = data_event_prob;
            event_vector_sum = nansum([event_vector_sum; ([(cosd(str2double(direction_name))), (sind(str2double(direction_name)))]*prob_)]);
            event_numTrials = event_numTrials + ...
                plot_data_compress.(['raster_data_' event_name]).(['train_data_logic_CS_' direction_name '_numTrial']);
        end
        event_tuning = atan2d(event_vector_sum(2), event_vector_sum(1));
        if (event_tuning<0); event_tuning=event_tuning+360; end
        plot_data_compress.CS_Tuning.([event_name '_tuning']) = event_tuning;
        overall_vector_sum = nansum([overall_vector_sum; event_vector_sum]);
        overall_prob = vertcat(overall_prob, plot_data_compress.CS_Tuning.([event_name '_prob']));
        overall_numTrials = vertcat(overall_numTrials, event_numTrials);
    end
    
    plot_data_compress.CS_Tuning.overall_prob = nanmean(overall_prob, 1);
    plot_data_compress.CS_Tuning.numTrials = nanmean(overall_numTrials, 1);
    % Using vector sum for CS_on
%     CS_on_index = round(ESN_Round(overall_tuning, 45) / 45); % CS_on_index is in 0-index format
%     overall_tuning = atan2d(overall_vector_sum(2), overall_vector_sum(1));
    % Using max prob for CS_on
    [~, CS_on_index] = max(plot_data_compress.CS_Tuning.overall_prob);
    CS_on_index = CS_on_index - 1; % CS_on_index should be in 0-index format
    overall_tuning = str2double(direction_names{CS_on_index+1});
    
    if (overall_tuning<0); overall_tuning=overall_tuning+360; end
    plot_data_compress.CS_Tuning.overall_tuning = overall_tuning;
    if (CS_on_index == 8); CS_on_index = 0; end
    plot_data_compress.CS_Tuning.overall_prob_tuning = mod((CS_on_index : 1 : CS_on_index+7), 8) + 1;
    
    %% Save plot_data_compress into ALL_PCELL_COMPRESSED_DATA
    ALL_PCELL_COMPRESSED_DATA(counter_pCell) = plot_data_compress; 
    fprintf(' --> Completed. \n')
end

end

%% function build_RASTER_DATA_ALL_PCELL_TUNED
function RASTER_DATA_ALL_PCELL_TUNED = build_RASTER_DATA_ALL_PCELL_TUNED(ALL_PCELL_COMPRESSED_DATA)
%% build RASTER_DATA_ALL_PCELL
fprintf(['Building RASTER_DATA_ALL_PCELL ' ' ... ']);
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED
clearvars RASTER_DATA_ALL_PCELL
for counter_ALL_PCELL = 1 : length(ALL_PCELL_COMPRESSED_DATA)
    field_names_ALL_PCELL = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL));
    field_names_ALL_PCELL = field_names_ALL_PCELL(contains(field_names_ALL_PCELL, 'raster_data_'));
    for counter_field_names_ALL_PCELL = 1 : 1 : length(field_names_ALL_PCELL)
        field_name_ALL_PCELL = field_names_ALL_PCELL{counter_field_names_ALL_PCELL};
        field_names_raster_data = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL));
        for counter_field_names_raster_data = 1 : 1 : length(field_names_raster_data)
            field_name_raster_data = field_names_raster_data{counter_field_names_raster_data};
            if contains(field_name_raster_data, '_sparse')
                continue;
            end
            RASTER_DATA_ALL_PCELL.(field_name_ALL_PCELL).(field_name_raster_data)(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).(field_name_raster_data);
        end
    end
end

train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
direction_names = {'000', '045', '090', '135', '180', '225', '270', '315'};
field_names_RASTER_DATA_ALL_PCELL = fieldnames(RASTER_DATA_ALL_PCELL);
for counter_field_names_RASTER_DATA_ALL_PCELL = 1 : 1 : length(field_names_RASTER_DATA_ALL_PCELL)
    field_name_RASTER_DATA_ALL_PCELL = field_names_RASTER_DATA_ALL_PCELL{counter_field_names_RASTER_DATA_ALL_PCELL};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        train_data_all = [];
        for counter_direction = 1 : 1 : length(direction_names)
            direction_name = direction_names{counter_direction};
            train_data_all = vertcat(train_data_all, RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name]));
        end
        RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all']) = train_data_all;
    end
    train_data_name = 'train_data_logic_CS';
    train_data_all_numTrial = [];
    for counter_direction = 1 : 1 : length(direction_names)
        direction_name = direction_names{counter_direction};
        train_data_all_numTrial = vertcat(train_data_all_numTrial, RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name '_numTrial']));
    end
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all_numTrial']) = train_data_all_numTrial;
    
    train_data_logic_SS_all = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all;
    train_data_logic_CS_all = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all;
    velocity_data_all       = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all;
    numTrial_all            = RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_numTrial;
    
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_mean = (numTrial_all' * train_data_logic_SS_all)./nansum(numTrial_all);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_mean = (numTrial_all' * train_data_logic_CS_all)./nansum(numTrial_all);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_mean       = (numTrial_all' * velocity_data_all)      ./nansum(numTrial_all);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_stdv = nanstd(train_data_logic_SS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_stdv = nanstd(train_data_logic_CS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_stdv       = nanstd(velocity_data_all,       numTrial_all, 1);
    
end

fprintf(' --> Completed. \n');

%% build RASTER_DATA_ALL_PCELL_TUNED
fprintf(['Building RASTER_DATA_ALL_PCELL_TUNED ' ' ... ']);
clearvars -except ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL RASTER_DATA_ALL_PCELL_TUNED
clearvars RASTER_DATA_ALL_PCELL_TUNED
train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
direction_names = {'000', '045', '090', '135', '180', '225', '270', '315'};

train_data_dirs  = {'_right', '_top', '_left', '_down'};
for counter_ALL_PCELL = 1 : length(ALL_PCELL_COMPRESSED_DATA)
    field_names_ALL_PCELL = fieldnames(ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL));
    field_names_ALL_PCELL = field_names_ALL_PCELL(contains(field_names_ALL_PCELL, 'raster_data_'));
    % Implement the CS-on direction with maximum CS probability
    overall_prob_tuning = ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).CS_Tuning.overall_prob_tuning;% 1:8;% 
    for counter_field_names_ALL_PCELL = 1 : 1 : length(field_names_ALL_PCELL)
        field_name_ALL_PCELL = field_names_ALL_PCELL{counter_field_names_ALL_PCELL};
        for counter_direction = 1 : 1 : length(direction_names)
            direction_name = direction_names{counter_direction};
            direction_name_tuned = direction_names{overall_prob_tuning(counter_direction)};
            for counter_train_data_names = 1 : 1 : length(train_data_names)
                train_data_name = train_data_names{counter_train_data_names};
                
                RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).([train_data_name '_' direction_name])(counter_ALL_PCELL, :) = ...
                    ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).([train_data_name '_' direction_name_tuned]);
            end
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['train_data_logic_CS_' direction_name '_numTrial'])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).(['train_data_logic_CS_' direction_name_tuned '_numTrial']);
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['numTrial_' direction_name])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).(['train_data_logic_CS_' direction_name_tuned '_numTrial']);
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['inds_span_' direction_name])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).('inds_span');
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['SS_firing_rate_' direction_name])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).('Neural_Properties_data').('SS_firing_rate');
            RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).(['CS_firing_rate_' direction_name])(counter_ALL_PCELL, :) = ...
                ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).('Neural_Properties_data').('CS_firing_rate');
        end
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_ALL_PCELL).('inds_span')(counter_ALL_PCELL, :) = ...
            ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL).(field_name_ALL_PCELL).('inds_span');
    end
end

train_data_names = {'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data', 'inds_span', 'numTrial', 'SS_firing_rate', 'CS_firing_rate'};
direction_names = {'000', '045', '090', '135', '180', '225', '270', '315'};
field_names_RASTER_DATA_ALL_PCELL = fieldnames(RASTER_DATA_ALL_PCELL_TUNED);
for counter_field_names_RASTER_DATA_ALL_PCELL = 1 : 1 : length(field_names_RASTER_DATA_ALL_PCELL)
    field_name_RASTER_DATA_ALL_PCELL = field_names_RASTER_DATA_ALL_PCELL{counter_field_names_RASTER_DATA_ALL_PCELL};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        
        train_data_all = [];
        for counter_direction = 1 : 1 : length(direction_names)
            direction_name = direction_names{counter_direction};
            train_data_all = vertcat(train_data_all, RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name]));
        end
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all']) = train_data_all;
        
        train_data_all = [];
        % CS_ON is defined as 315, 000, and 045
        for counter_direction = [8 1 2]
            direction_name = direction_names{counter_direction};
            train_data_all = vertcat(train_data_all, RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name]));
        end
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_ON']) = train_data_all;
        
        train_data_all = [];
        % CS_OFF is defined as 135, 180, and 225
        for counter_direction = [4 5 6]
            direction_name = direction_names{counter_direction};
            train_data_all = vertcat(train_data_all, RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name]));
        end
        RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_OFF']) = train_data_all;
    end
    train_data_name = 'train_data_logic_CS';
    train_data_all_numTrial = [];
    for counter_direction = 1 : 1 : length(direction_names)
        direction_name = direction_names{counter_direction};
        train_data_all_numTrial = vertcat(train_data_all_numTrial, RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name '_numTrial']));
    end
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_all_numTrial']) = train_data_all_numTrial;
    
    train_data_all_numTrial = [];
    % CS_ON is defined as 315, 000, and 045
    for counter_direction = [8 1 2]
        direction_name = direction_names{counter_direction};
        train_data_all_numTrial = vertcat(train_data_all_numTrial, RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name '_numTrial']));
    end
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_ON_numTrial']) = train_data_all_numTrial;
    
    train_data_all_numTrial = [];
    % CS_OFF is defined as 135, 180, and 225
    for counter_direction = [4 5 6]
        direction_name = direction_names{counter_direction};
        train_data_all_numTrial = vertcat(train_data_all_numTrial, RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([train_data_name '_' direction_name '_numTrial']));
    end
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).([    train_data_name '_OFF_numTrial']) = train_data_all_numTrial;
    
    train_data_logic_SS_all = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all;
    train_data_logic_CS_all = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all;
    velocity_data_all       = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all;
    numTrial_all            = RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_numTrial;
    
    train_data_logic_SS_all(isnan(train_data_logic_SS_all)) = 0;
    train_data_logic_CS_all(isnan(train_data_logic_CS_all)) = 0;
    velocity_data_all(isnan(velocity_data_all)) = 0;
    numTrial_all(isnan(numTrial_all)) = 0;
    
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_mean = (numTrial_all' * train_data_logic_SS_all)./nansum(numTrial_all);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_mean = (numTrial_all' * train_data_logic_CS_all)./nansum(numTrial_all);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_mean       = (numTrial_all' * velocity_data_all)      ./nansum(numTrial_all);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_SS_all_stdv = nanstd(train_data_logic_SS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).train_data_logic_CS_all_stdv = nanstd(train_data_logic_CS_all, numTrial_all, 1);
    RASTER_DATA_ALL_PCELL_TUNED.(field_name_RASTER_DATA_ALL_PCELL).velocity_data_all_stdv       = nanstd(velocity_data_all,       numTrial_all, 1);
    
end

fprintf(' --> Completed. \n');

end

%% function plot_RASTER_DATA_ALL_PCELL_TUNED
function plot_RASTER_DATA_ALL_PCELL_TUNED(path_data_monkey_sorted, ALL_PCELL_COMPRESSED_DATA, RASTER_DATA_ALL_PCELL_TUNED)
%% Plot RASTER_DATA_ALL_PCELL_TUNED
fprintf(['Ploting ' ' ... ']);
close all
range_SS = [20 120];
range_SS_population = [-40 60];
range_CS = [0 4];
range_vm = [0 600];
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
clearvars hFig hAx

train_data_dirs  = {'_000', '_045', '_090', '_135', '_180', '_225', '_270', '_315', '_all', '_ON', '_OFF'};
train_data_names = {'train_data_logic_SS', 'train_data_logic_SS', 'train_data_logic_CS', 'velocity_data'};
field_names = {'raster_data_cue_present', 'raster_data_primSac_onset', 'raster_data_primSac_offset', 'raster_data_corrSac_onset' };
yLabels = {'Change in Population (spk/s)', 'SS Firing (spk/s)', 'CS Firing (spk/s)', 'Eye velocity (deg/s)'};
xLabels = {'Cue Pres (ms)', 'Prim Sac On (ms)', 'Prim Sac Off (ms)', 'Corr Sac On (ms)'};
sgTitles = {'CS-000', 'CS-045', 'CS-090', 'CS-135', 'CS-180', 'CS-225', 'CS-270', 'CS-315', 'All Direction', 'CS-ON', 'CS-OFF'};
for counter_train_data_dirs = 1 : 1 : length(train_data_dirs)
    hFig(counter_train_data_dirs) = figure(counter_train_data_dirs);
    clf(hFig(counter_train_data_dirs));
    train_data_dir = train_data_dirs{counter_train_data_dirs};
    for counter_train_data_names = 1 : 1 : length(train_data_names)
        train_data_name = train_data_names{counter_train_data_names};
        for counter_field_names = 1 : 1 : length(field_names)
            field_name = field_names{counter_field_names};
            inds_span = nanmean(RASTER_DATA_ALL_PCELL_TUNED.(field_name).inds_span, 1);
            if counter_train_data_names == 1
                numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
                SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
                train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
                train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);
                
                num_perm = 5000;
                num_samples = num_pCells ; % 50; % 
                length_span = size(train_data_logic_SS, 2);
                train_data_logic_SS_perm = nan(num_perm, length_span);
                for counter_perm = 1 : num_perm
                    inds_perm = randi(num_pCells,[num_samples 1]);
                    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
                    train_data_logic_SS_(isnan(train_data_logic_SS_)) = 0;
                    numTrial_ = numTrial(inds_perm, :);
                    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ nansum(numTrial_) .* 1000;
                    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
                end
                variable_data_mean = ESN_smooth( nanmean(train_data_logic_SS_perm, 1) );
                variable_data_stdv = ESN_smooth( nanstd(train_data_logic_SS_perm, 0, 1) ) ./ sqrt(num_samples) * 3; % 95% CI
            end
            if counter_train_data_names ~= 1
                variable_data = RASTER_DATA_ALL_PCELL_TUNED.(field_name).([train_data_name train_data_dir]);
                variable_data_mean = nanmean(variable_data, 1);
                variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells);
                if(counter_train_data_names ~= 4)
                    variable_data_mean = ESN_smooth(variable_data_mean * 1000);
                    variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
                end
            end
            
            subplot_num = 4*(counter_train_data_names-1)+counter_field_names;
            hAx{counter_train_data_dirs}(subplot_num) = subplot(4, 4, subplot_num);
            hAx_ = hAx{counter_train_data_dirs}(subplot_num);
            hold(hAx_, 'on')
            plot(inds_span, variable_data_mean+variable_data_stdv, '-k', 'LineWidth', 0.5)
            plot(inds_span, variable_data_mean-variable_data_stdv, '-k', 'LineWidth', 0.5)
            plot(inds_span, variable_data_mean, '-k', 'LineWidth', 1)
            if(counter_field_names == 1)
                ylabel(hAx_, yLabels{counter_train_data_names});
            end
            if(counter_train_data_names == 4)
                xlabel(hAx_, xLabels{counter_field_names});
            end
            if(counter_train_data_names == 1)
                ylim(range_SS_population);
            elseif(counter_train_data_names == 2)
                ylim(range_SS);
            elseif(counter_train_data_names == 3)
                ylim(range_CS);
            elseif(counter_train_data_names == 4)
                ylim(range_vm);
            end
        end
    end
    ESN_Beautify_Plot(hFig(counter_train_data_dirs), [11, 8.5])
    sgtitle(hFig(counter_train_data_dirs), sgTitles{counter_train_data_dirs}, 'Interpreter', 'none');
end

fprintf(' --> Completed. \n');

%% plot CS_Tune
fprintf(['Ploting ' ' ... ']);
hFig(length(train_data_dirs)+1) = figure(length(train_data_dirs)+1);
clf(hFig(length(train_data_dirs)+1));
cue_present_prob_ALL_PCELL    = nan(length(ALL_PCELL_COMPRESSED_DATA), 8);
primSac_onset_prob_ALL_PCELL  = nan(length(ALL_PCELL_COMPRESSED_DATA), 8);
primSac_offset_prob_ALL_PCELL = nan(length(ALL_PCELL_COMPRESSED_DATA), 8);
corrSac_onset_prob_ALL_PCELL  = nan(length(ALL_PCELL_COMPRESSED_DATA), 8);
numTrials_ALL_PCELL = nan(length(ALL_PCELL_COMPRESSED_DATA), 1);
for counter_ALL_PCELL_COMPRESSED_DATA = 1 : 1 : length(ALL_PCELL_COMPRESSED_DATA)
    % Implement the CS-on direction with maximum CS probability
    overall_prob_tuning = ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.overall_prob_tuning;% 1:8;% 
    cue_present_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.cue_present_prob( overall_prob_tuning );
    primSac_onset_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.primSac_onset_prob( overall_prob_tuning );
    primSac_offset_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.primSac_offset_prob( overall_prob_tuning );
    corrSac_onset_prob_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ...
        ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.corrSac_onset_prob(overall_prob_tuning);
    numTrials_ALL_PCELL(counter_ALL_PCELL_COMPRESSED_DATA, :) = ALL_PCELL_COMPRESSED_DATA(counter_ALL_PCELL_COMPRESSED_DATA).CS_Tuning.numTrials;
end

num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
cue_present_prob_TUNED    = nan(num_pCells, 8);
primSac_onset_prob_TUNED  = nan(num_pCells, 8);
primSac_offset_prob_TUNED = nan(num_pCells, 8);
corrSac_onset_prob_TUNED  = nan(num_pCells, 8);
for counter_pCell = 1 : num_pCells
    cue_present_prob_TUNED(counter_pCell, :) = ...
        numTrials_ALL_PCELL(counter_pCell, :)' * cue_present_prob_ALL_PCELL(counter_pCell, :) ./ nansum(numTrials_ALL_PCELL(counter_pCell, :));
    primSac_onset_prob_TUNED(counter_pCell, :) = ...
        numTrials_ALL_PCELL(counter_pCell, :)' * primSac_onset_prob_ALL_PCELL(counter_pCell, :) ./ nansum(numTrials_ALL_PCELL(counter_pCell, :));
    primSac_offset_prob_TUNED(counter_pCell, :) = ...
        numTrials_ALL_PCELL(counter_pCell, :)' * primSac_offset_prob_ALL_PCELL(counter_pCell, :) ./ nansum(numTrials_ALL_PCELL(counter_pCell, :));
    corrSac_onset_prob_TUNED(counter_pCell, :) = ...
        numTrials_ALL_PCELL(counter_pCell, :)' * corrSac_onset_prob_ALL_PCELL(counter_pCell, :) ./ nansum(numTrials_ALL_PCELL(counter_pCell, :));
end

overall_prob_TUNED = (cue_present_prob_TUNED + ...
    primSac_onset_prob_TUNED + ...
    primSac_offset_prob_TUNED + ...
    corrSac_onset_prob_TUNED) ./ 4;

overall_prob_TUNED_mean = nanmean(overall_prob_TUNED, 1);
overall_prob_TUNED_stdv = nanstd(overall_prob_TUNED, 0, 1) ./ sqrt(num_pCells);
overall_prob_TUNED_stdv_plus = overall_prob_TUNED_mean + overall_prob_TUNED_stdv;
overall_prob_TUNED_stdv_minus = overall_prob_TUNED_mean - overall_prob_TUNED_stdv;

plot_data_amp_mean = [overall_prob_TUNED_mean, overall_prob_TUNED_mean(1), nan]';
plot_data_deg_mean = [0, 45, 90, 135, 180, 225, 270, 315, 0, nan]';
plot_data_x_axis_mean = plot_data_amp_mean .* cosd(plot_data_deg_mean);
plot_data_y_axis_mean = plot_data_amp_mean .* sind(plot_data_deg_mean);

plot_data_amp_stdv_p = [overall_prob_TUNED_stdv_plus, overall_prob_TUNED_stdv_plus(1), nan]';
plot_data_deg_stdv_p = [0, 45, 90, 135, 180, 225, 270, 315, 0, nan]';
plot_data_x_axis_stdv_p = plot_data_amp_stdv_p .* cosd(plot_data_deg_stdv_p);
plot_data_y_axis_stdv_p = plot_data_amp_stdv_p .* sind(plot_data_deg_stdv_p);

plot_data_amp_stdv_m = [overall_prob_TUNED_stdv_minus, overall_prob_TUNED_stdv_minus(1), nan]';
plot_data_deg_stdv_m = [0, 45, 90, 135, 180, 225, 270, 315, 0, nan]';
plot_data_x_axis_stdv_m = plot_data_amp_stdv_m .* cosd(plot_data_deg_stdv_m);
plot_data_y_axis_stdv_m = plot_data_amp_stdv_m .* sind(plot_data_deg_stdv_m);

hold on
plot([0.4, -0.4, nan, 0, -0, nan,]', [0, 0, nan, 0.4, -0.4, nan], 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
plot(cosd(0:5:360)*0.05, sind(0:5:360)*0.05, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.15, sind(0:5:360)*0.15, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.25, sind(0:5:360)*0.25, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.8 0.8 0.8])

plot_data_x_axis_mean   = plot_data_x_axis_mean(  ~isnan(plot_data_x_axis_mean)); 
plot_data_y_axis_mean   = plot_data_y_axis_mean(  ~isnan(plot_data_y_axis_mean));
plot_data_x_axis_stdv_p = plot_data_x_axis_stdv_p(~isnan(plot_data_x_axis_stdv_p)); 
plot_data_y_axis_stdv_p = plot_data_y_axis_stdv_p(~isnan(plot_data_y_axis_stdv_p));
plot_data_x_axis_stdv_m = plot_data_x_axis_stdv_m(~isnan(plot_data_x_axis_stdv_m)); 
plot_data_y_axis_stdv_m = plot_data_y_axis_stdv_m(~isnan(plot_data_y_axis_stdv_m));

plot(plot_data_x_axis_mean,   plot_data_y_axis_mean, '-k', 'LineWidth', 1)
plot(plot_data_x_axis_stdv_p, plot_data_y_axis_stdv_p, '-k', 'LineWidth', 0.5)
plot(plot_data_x_axis_stdv_m, plot_data_y_axis_stdv_m, '-k', 'LineWidth', 0.5)

axis equal
xlim([-0.3 0.3])
ylim([-0.3 0.3])
xlabel('CS probability');
ylabel('CS probability');
set(gca, 'XTick', -0.3:0.1:0.3, 'YTick', -0.3:0.1:0.3)

ESN_Beautify_Plot(hFig(length(train_data_dirs)+1), [4, 3])
title('CS Tuning', 'Interpreter', 'none');

hFig(length(train_data_dirs)+2) = figure(length(train_data_dirs)+2);
clf(hFig(length(train_data_dirs)+2));
hold on
plot_order_ = [7 8 1 2 3 4 5 6];
plot(overall_prob_TUNED_stdv_plus(plot_order_), '-k', 'LineWidth', 0.5)
plot(overall_prob_TUNED_stdv_minus(plot_order_), '-k', 'LineWidth', 0.5)
plot(overall_prob_TUNED_mean(plot_order_), '-k', 'LineWidth', 1)
ylabel('CS probability');
xlabel('Direction')
set(gca, 'XTick', 1:1:8, 'XTickLabel', {'-90','-45','ON','45','90','135','180','225'})

avg_prob_TUNED_mean = nanmean(nanmean(overall_prob_TUNED, 2));
avg_prob_TUNED_stdv = nanstd(nanmean(overall_prob_TUNED, 2), 0, 1) ./ sqrt(num_pCells);
avg_prob_TUNED_mean = repmat(avg_prob_TUNED_mean, 1, size(overall_prob_TUNED,2));
avg_prob_TUNED_stdv = repmat(avg_prob_TUNED_stdv, 1, size(overall_prob_TUNED,2));
avg_prob_TUNED_stdv_plus = avg_prob_TUNED_mean + avg_prob_TUNED_stdv;
avg_prob_TUNED_stdv_minus = avg_prob_TUNED_mean - avg_prob_TUNED_stdv;

plot_order_ = [7 8 1 2 3 4 5 6];
plot(avg_prob_TUNED_stdv_plus(plot_order_), '-k', 'LineWidth', 0.5)
plot(avg_prob_TUNED_stdv_minus(plot_order_), '-k', 'LineWidth', 0.5)
plot(avg_prob_TUNED_mean(plot_order_), '-k', 'LineWidth', 1)
ylim([0.10 0.30])

ESN_Beautify_Plot(hFig(length(train_data_dirs)+2), [4, 3])
title('CS Tuning', 'Interpreter', 'none');

fprintf(' --> Completed. \n');

%% Save combined file
fprintf(['Saving plots ' ' ... ']);
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
save_file_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
save_file_path = path_data_monkey_sorted;
if ~exist(save_file_path, 'dir')
    mkdir(save_file_path);
end
saveas(hFig(1), [save_file_path filesep save_file_name '_CS_000'], 'pdf');
saveas(hFig(1), [save_file_path filesep save_file_name '_CS_000'], 'png');
saveas(hFig(2), [save_file_path filesep save_file_name '_CS_045'], 'pdf');
saveas(hFig(2), [save_file_path filesep save_file_name '_CS_045'], 'png');
saveas(hFig(3), [save_file_path filesep save_file_name '_CS_090'], 'pdf');
saveas(hFig(3), [save_file_path filesep save_file_name '_CS_090'], 'png');
saveas(hFig(4), [save_file_path filesep save_file_name '_CS_135'], 'pdf');
saveas(hFig(4), [save_file_path filesep save_file_name '_CS_135'], 'png');
saveas(hFig(5), [save_file_path filesep save_file_name '_CS_180'], 'pdf');
saveas(hFig(5), [save_file_path filesep save_file_name '_CS_180'], 'png');
saveas(hFig(6), [save_file_path filesep save_file_name '_CS_225'], 'pdf');
saveas(hFig(6), [save_file_path filesep save_file_name '_CS_225'], 'png');
saveas(hFig(7), [save_file_path filesep save_file_name '_CS_270'], 'pdf');
saveas(hFig(7), [save_file_path filesep save_file_name '_CS_270'], 'png');
saveas(hFig(8), [save_file_path filesep save_file_name '_CS_315'], 'pdf');
saveas(hFig(8), [save_file_path filesep save_file_name '_CS_315'], 'png');
saveas(hFig(9), [save_file_path filesep save_file_name '_CS_all'], 'pdf');
saveas(hFig(9), [save_file_path filesep save_file_name '_CS_all'], 'png');
saveas(hFig(10), [save_file_path filesep save_file_name '_CS_ON' ], 'pdf');
saveas(hFig(10), [save_file_path filesep save_file_name '_CS_ON' ], 'png');
saveas(hFig(11), [save_file_path filesep save_file_name '_CS_OFF'], 'pdf');
saveas(hFig(11), [save_file_path filesep save_file_name '_CS_OFF'], 'png');
saveas(hFig(12), [save_file_path filesep save_file_name '_CS_tun'], 'pdf');
saveas(hFig(12), [save_file_path filesep save_file_name '_CS_tun'], 'png');
saveas(hFig(13), [save_file_path filesep save_file_name '_CS_tun_linear'], 'pdf');
saveas(hFig(13), [save_file_path filesep save_file_name '_CS_tun_linear'], 'png');

close all

fprintf(' --> Completed. \n');

end

%% function ESN_smooth
function smooth_data_ = ESN_smooth(data_)
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
smooth_data_ = smooth(data_, 21, 'sgolay', 2);
end

%% function RE-RUN ESN_plot_neural_modulation_v4
function re_run_ESN_plot_neural_modulation_v4(pCell_list, path_data_monkey_sorted)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name(1:2);
        month_ = file_name(3:4);
        day_ = file_name(5:6);
        hour_ = file_name(8:9);
        minute_ = file_name(10:11);
        second_ = file_name(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        %% RE-RUN ESN_plot_neural_modulation_v4
        params.auto = true;
        params.file_path = file_path;
        params.file_name = file_name;
        ESN_plot_neural_modulation_v4(1, params);
    end
end

end

%% function RE-RUN ESN_monkey_behavior_calibrate
function re_run_ESN_monkey_behavior_calibrate(pCell_list, path_data_monkey_sorted)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    file_path_list = cell(0,1);
    file_name_list = cell(0,1);
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name(1:2);
        month_ = file_name(3:4);
        day_ = file_name(5:6);
        hour_ = file_name(8:9);
        minute_ = file_name(10:11);
        second_ = file_name(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        file_name_ANALYZED = [file_name(1:13) '_ANALYZED.mat'];
        file_path_list{counter_recording,1} = file_path;
        file_name_list{counter_recording,1} = file_name_ANALYZED;
    end
    %% RE-RUN ESN_monkey_behavior_calibrate
    ESN_monkey_behavior_calibrate(file_path_list, file_name_list);
end

end

%% function plot_SS_CS_firing_histogram
function plot_SS_CS_firing_histogram()
clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
numTrials_pCells = nan(num_pCells, 1);
SS_firing_pCells = nan(num_pCells, 1);
CS_firing_pCells = nan(num_pCells, 1);
for counter_pCell = 1 : 1 : num_pCells
    numTrials_pCells(counter_pCell, 1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).CS_Tuning.numTrials;
    SS_firing_pCells(counter_pCell, 1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties_data.SS_firing_rate;
    CS_firing_pCells(counter_pCell, 1) = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties_data.CS_firing_rate;
end

SS_firing_pCells_mean = nanmean(SS_firing_pCells);
SS_firing_pCells_stdv = nanstd( SS_firing_pCells);
SS_firing_pCells_sem  = nanstd( SS_firing_pCells)./sqrt(num_pCells);
CS_firing_pCells_mean = nanmean(CS_firing_pCells);
CS_firing_pCells_stdv = nanstd( CS_firing_pCells);
CS_firing_pCells_sem  = nanstd( CS_firing_pCells)./sqrt(num_pCells);

fprintf(['numTrials ', 'mean: ', num2str(mean(numTrials_pCells)), ',\t SEM: ', num2str(std(numTrials_pCells)./sqrt(num_pCells)), ',\t std: ', num2str(std(numTrials_pCells)), '\n'])
fprintf(['SS_firing ', 'mean: ', num2str(mean(SS_firing_pCells)), ',\t SEM: ', num2str(std(SS_firing_pCells)./sqrt(num_pCells)), ',\t std: ', num2str(std(SS_firing_pCells)), '\n'])
fprintf(['CS_firing ', 'mean: ', num2str(mean(CS_firing_pCells)), ',\t SEM: ', num2str(std(CS_firing_pCells)./sqrt(num_pCells)), ',\t std: ', num2str(std(CS_firing_pCells)), '\n'])

hFig = figure(101);
clf(hFig)
subplot(1,2,1)
hold on
SS_firing_edges = 5:10:135;
histogram(SS_firing_pCells, SS_firing_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'b')
histogram(SS_firing_pCells, SS_firing_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([(SS_firing_pCells_mean+SS_firing_pCells_stdv), (SS_firing_pCells_mean+SS_firing_pCells_stdv)], yl_, '--b', 'LineWidth', 0.5);
plot([(SS_firing_pCells_mean-SS_firing_pCells_stdv), (SS_firing_pCells_mean-SS_firing_pCells_stdv)], yl_, '--b', 'LineWidth', 0.5);
plot([(SS_firing_pCells_mean+SS_firing_pCells_sem), (SS_firing_pCells_mean+SS_firing_pCells_sem)], yl_, '--b', 'LineWidth', 1);
plot([(SS_firing_pCells_mean-SS_firing_pCells_sem), (SS_firing_pCells_mean-SS_firing_pCells_sem)], yl_, '--b', 'LineWidth', 1);
plot([(SS_firing_pCells_mean), (SS_firing_pCells_mean)], yl_, '-b', 'LineWidth', 1.5);
xlabel('SS Firing Rate (Hz)')
ylabel('Count (#)')
title(['mean= ' num2str(SS_firing_pCells_mean,3)])

subplot(1,2,2)
hold on
CS_firing_edges = 0.25:0.1:1.55;
histogram(CS_firing_pCells, CS_firing_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(CS_firing_pCells, CS_firing_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([(CS_firing_pCells_mean+CS_firing_pCells_stdv), (CS_firing_pCells_mean+CS_firing_pCells_stdv)], yl_, '--r', 'LineWidth', 0.5);
plot([(CS_firing_pCells_mean-CS_firing_pCells_stdv), (CS_firing_pCells_mean-CS_firing_pCells_stdv)], yl_, '--r', 'LineWidth', 0.5);
plot([(CS_firing_pCells_mean+CS_firing_pCells_sem),  (CS_firing_pCells_mean+CS_firing_pCells_sem)],  yl_, '--r', 'LineWidth', 1);
plot([(CS_firing_pCells_mean-CS_firing_pCells_sem),  (CS_firing_pCells_mean-CS_firing_pCells_sem)],  yl_, '--r', 'LineWidth', 1);
plot([(CS_firing_pCells_mean),                       (CS_firing_pCells_mean)], yl_, '-r', 'LineWidth', 1.5);
xlabel('CS Firing Rate (Hz)')
ylabel('Count (#)')
title(['mean= ' num2str(CS_firing_pCells_mean,2)])
ESN_Beautify_Plot(hFig, [4, 2], 8)

end

%% function plot_waveforms
function plot_waveforms()
clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
numTrials_pCells = nan(num_pCells, 1);
SS_waveform_pCells = nan(num_pCells, 180);
CS_waveform_pCells = nan(num_pCells, 180);
SS_xprob_pCells = nan(num_pCells, 100);
CS_xprob_pCells = nan(num_pCells, 100);
SS_firing_pCells = nan(num_pCells, 1);
for counter_pCell = 1 : 1 : num_pCells
    Neural_Properties_data = ALL_PCELL_COMPRESSED_DATA(counter_pCell).Neural_Properties_data;
    SS_waveform_pCells(counter_pCell, :) = Neural_Properties_data.SS_waveform;
    CS_waveform_pCells(counter_pCell, :) = Neural_Properties_data.CS_waveform;
    SS_xprob_pCells(counter_pCell, :) = Neural_Properties_data.Corr_data.SS_SSxSS_AUTO;
    CS_xprob_pCells(counter_pCell, :) = Neural_Properties_data.Corr_data.CS_CSxSS_AUTO;
    SS_firing_pCells(counter_pCell, 1) = Neural_Properties_data.SS_firing_rate;
end
SS_waveform_pCells_max = max([max(SS_waveform_pCells, [], 2), abs(min(SS_waveform_pCells, 2))], [], 2);
SS_waveform_pCells_norm = SS_waveform_pCells ./ SS_waveform_pCells_max;
CS_waveform_pCells_norm = CS_waveform_pCells ./ SS_waveform_pCells_max;
SS_xprob_pCells_norm = SS_xprob_pCells ./ (SS_firing_pCells./1000);
SS_xprob_pCells_norm(:,50) = 0;
CS_xprob_pCells_norm = CS_xprob_pCells ./ (SS_firing_pCells./1000);

SS_waveform_mean = nanmean(SS_waveform_pCells_norm);
SS_waveform_stdv = nanstd( SS_waveform_pCells_norm);%./sqrt(num_pCells);
SS_waveform_stdv_p = SS_waveform_mean + SS_waveform_stdv;
SS_waveform_stdv_m = SS_waveform_mean - SS_waveform_stdv;

CS_waveform_mean = nanmean(CS_waveform_pCells_norm);
CS_waveform_stdv = nanstd( CS_waveform_pCells_norm);%./sqrt(num_pCells);
CS_waveform_stdv_p = CS_waveform_mean + CS_waveform_stdv;
CS_waveform_stdv_m = CS_waveform_mean - CS_waveform_stdv;

SS_xprob_mean = nanmean(SS_xprob_pCells_norm);
SS_xprob_stdv = nanstd( SS_xprob_pCells_norm);%./sqrt(num_pCells);
SS_xprob_stdv_p = SS_xprob_mean + SS_xprob_stdv;
SS_xprob_stdv_m = SS_xprob_mean - SS_xprob_stdv;

CS_xprob_mean = nanmean(CS_xprob_pCells_norm);
CS_xprob_stdv = nanstd( CS_xprob_pCells_norm);%./sqrt(num_pCells);
CS_xprob_stdv_p = CS_xprob_mean + CS_xprob_stdv;
CS_xprob_stdv_m = CS_xprob_mean - CS_xprob_stdv;

time_waveform = (-59:1:120)./30;
time_xprob = (-49:1:50);

hFig = figure(105);
clf(hFig)
subplot(1,2,1)
hold on
% plot(time_waveform, SS_waveform_pCells_norm, '-k', 'LineWidth', 0.5)
plot(time_waveform, SS_waveform_stdv_m, '-b', 'LineWidth', 0.5)
plot(time_waveform, SS_waveform_stdv_p, '-b', 'LineWidth', 0.5)
plot(time_waveform, SS_waveform_mean, '-b', 'LineWidth', 1.5)
% ylabel('SS waveform')
% xlabel('Time (ms)')

% subplot(2,2,2)
% hold on
% plot(time_waveform, CS_waveform_pCells_norm, '-k', 'LineWidth', 0.5)
plot(time_waveform, CS_waveform_stdv_m, '-r', 'LineWidth', 0.5)
plot(time_waveform, CS_waveform_stdv_p, '-r', 'LineWidth', 0.5)
plot(time_waveform, CS_waveform_mean, '-r', 'LineWidth', 1.5)
ylabel('waveform')
xlabel('Time (ms)')


subplot(1,2,2)
hold on
% plot(time_xprob, SS_xprob_pCells_norm, '-k', 'LineWidth', 0.5)
plot(time_xprob, SS_xprob_stdv_p, '-b', 'LineWidth', 0.5)
plot(time_xprob, SS_xprob_stdv_m, '-b', 'LineWidth', 0.5)
plot(time_xprob, SS_xprob_mean, '-b', 'LineWidth', 1.5)
% ylabel('SS prob')
% xlabel('Time (ms)')

% subplot(2,2,4)
% hold on
% plot(time_xprob, CS_xprob_pCells_norm, '-k', 'LineWidth', 0.5)
plot(time_xprob, CS_xprob_stdv_p, '-r', 'LineWidth', 0.5)
plot(time_xprob, CS_xprob_stdv_m, '-r', 'LineWidth', 0.5)
plot(time_xprob, CS_xprob_mean, '-r', 'LineWidth', 1.5)
ylabel('prob')
xlabel('Time (ms)')

ESN_Beautify_Plot(hFig, [4, 2], 8)

%%
hFig = figure(106);
clf(hFig)
% time_xprob = (-49:1:50);
CS_xprob_pCells_norm = CS_xprob_pCells ./ nanmean(CS_xprob_pCells(:,20:50), 2);
% CS_xprob_pCells_norm = smoothdata(CS_xprob_pCells_norm, 2, 'sgolay', 5, 'omitnan', 'Degree', 2);
[~,idx] = max(CS_xprob_pCells_norm(:,56:end)>0.63, [], 2);
suppression_ = idx+5;
hold on
edges_ = 5.5: 1 : 25.5;
histogram(suppression_, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(suppression_, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 2)
xline(nanmean(suppression_), 'r', 'linewidth', 1)
ylabel('Count')
xlabel('CS suppression (ms)')
ESN_Beautify_Plot(hFig, [2, 1.5], 8)

end

%% function generate_pCell_id_for_excel
function generate_pCell_id_for_excel
%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
pCell_id_for_excel = cell(0,1);
for counter_pCell = 1 : 1 : num_pCells
    id_ = ALL_PCELL_COMPRESSED_DATA(counter_pCell).id;
    if length(id_) == 1
        name_ = id_{1};
    else
        name_ = [id_{1} '_combined_' num2str(length(id_))];
    end
    pCell_id_for_excel(counter_pCell, 1) = {name_};
end

end

%% function plot_DCN
function plot_DCN
clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
inds_span = nanmean(RASTER_DATA_ALL_PCELL_TUNED.('raster_data_primSac_onset').inds_span, 1);
%
% SS_prim_000 = RASTER_DATA_ALL_PCELL_TUNED.('raster_data_primSac_onset').(['train_data_logic_SS' '_000']);
% SS_prim_180 = RASTER_DATA_ALL_PCELL_TUNED.('raster_data_primSac_onset').(['train_data_logic_SS' '_180']);
% SS_prim_DCN = SS_prim_180 - SS_prim_000;

field_name = 'raster_data_primSac_onset';
train_data_dir = '_180';
numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);
num_perm = 5000;
num_samples = num_pCells ; % 50; %
length_span = size(train_data_logic_SS, 2);
train_data_logic_SS_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
    train_data_logic_SS_(isnan(train_data_logic_SS_)) = 0;
    numTrial_ = numTrial(inds_perm, :);
    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ nansum(numTrial_);
    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
end
SS_prim_180 = train_data_logic_SS_perm;

field_name = 'raster_data_primSac_onset';
train_data_dir = '_000';
numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);
num_perm = 5000;
num_samples = num_pCells ; % 50; %
length_span = size(train_data_logic_SS, 2);
train_data_logic_SS_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
    train_data_logic_SS_(isnan(train_data_logic_SS_)) = 0;
    numTrial_ = numTrial(inds_perm, :);
    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ nansum(numTrial_);
    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
end
SS_prim_000 = train_data_logic_SS_perm;

SS_prim_DCN = SS_prim_180 - SS_prim_000;

%
% SS_corr_000 = RASTER_DATA_ALL_PCELL_TUNED.('raster_data_corrSac_onset').(['train_data_logic_SS' '_000']);
% SS_corr_180 = RASTER_DATA_ALL_PCELL_TUNED.('raster_data_corrSac_onset').(['train_data_logic_SS' '_180']);
% SS_corr_DCN = SS_corr_180 - SS_corr_000;

field_name = 'raster_data_corrSac_onset';
train_data_dir = '_180';
numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);
num_perm = 5000;
num_samples = num_pCells ; % 50; %
length_span = size(train_data_logic_SS, 2);
train_data_logic_SS_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
    train_data_logic_SS_(isnan(train_data_logic_SS_)) = 0;
    numTrial_ = numTrial(inds_perm, :);
    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ nansum(numTrial_);
    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
end
SS_corr_180 = train_data_logic_SS_perm;

field_name = 'raster_data_corrSac_onset';
train_data_dir = '_000';
numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);
num_perm = 5000;
num_samples = num_pCells ; % 50; %
length_span = size(train_data_logic_SS, 2);
train_data_logic_SS_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
    train_data_logic_SS_(isnan(train_data_logic_SS_)) = 0;
    numTrial_ = numTrial(inds_perm, :);
    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ nansum(numTrial_);
    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
end
SS_corr_000 = train_data_logic_SS_perm;

SS_corr_DCN = SS_corr_180 - SS_corr_000;
%
range_SS = [-15 25];
range_SS_DCN = [-15 25];
range_idx = 200:400;
hFig = figure(102);
clf(hFig);
subplot(2,3,1)
hold on
variable_data = SS_prim_000;
variable_data_mean = nanmean(variable_data, 1);
variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells) * 3; % 95% CI
variable_data_mean = ESN_smooth(variable_data_mean * 1000);
variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
plot(inds_span(range_idx), variable_data_mean(range_idx)+variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx)-variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx), '-k', 'LineWidth', 1)
ylim(range_SS)
title('SS_prim_000', 'Interpreter', 'none')

subplot(2,3,2)
hold on
variable_data = SS_prim_180;
variable_data_mean = nanmean(variable_data, 1);
variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells) * 3; % 95% CI
variable_data_mean = ESN_smooth(variable_data_mean * 1000);
variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
plot(inds_span(range_idx), variable_data_mean(range_idx)+variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx)-variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx), '-k', 'LineWidth', 1)
ylim(range_SS)
title('SS_prim_180', 'Interpreter', 'none')

subplot(2,3,3)
hold on
variable_data = SS_prim_DCN;
variable_data_mean = nanmean(variable_data, 1);
variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells) * 3; % 95% CI
variable_data_mean = ESN_smooth(variable_data_mean * 1000);
variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
plot(inds_span(range_idx), variable_data_mean(range_idx)+variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx)-variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx), '-k', 'LineWidth', 1)
ylim(range_SS_DCN)
title('SS_prim_DCN', 'Interpreter', 'none')

subplot(2,3,4)
hold on
variable_data = SS_corr_000;
variable_data_mean = nanmean(variable_data, 1);
variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells) * 3; % 95% CI
variable_data_mean = ESN_smooth(variable_data_mean * 1000);
variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
plot(inds_span(range_idx), variable_data_mean(range_idx)+variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx)-variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx), '-k', 'LineWidth', 1)
ylim(range_SS)
title('SS_corr_000', 'Interpreter', 'none')

subplot(2,3,5)
hold on
variable_data = SS_corr_180;
variable_data_mean = nanmean(variable_data, 1);
variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells) * 3; % 95% CI
variable_data_mean = ESN_smooth(variable_data_mean * 1000);
variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
plot(inds_span(range_idx), variable_data_mean(range_idx)+variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx)-variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx), '-k', 'LineWidth', 1)
ylim(range_SS)
title('SS_corr_180', 'Interpreter', 'none')

subplot(2,3,6)
hold on
variable_data = SS_corr_DCN;
variable_data_mean = nanmean(variable_data, 1);
variable_data_stdv = nanstd(variable_data, 0, 1) ./ sqrt(num_pCells) * 3; % 95% CI
variable_data_mean = ESN_smooth(variable_data_mean * 1000);
variable_data_stdv = ESN_smooth(variable_data_stdv * 1000);
plot(inds_span(range_idx), variable_data_mean(range_idx)+variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx)-variable_data_stdv(range_idx), '-k', 'LineWidth', 0.5)
plot(inds_span(range_idx), variable_data_mean(range_idx), '-k', 'LineWidth', 1)
ylim(range_SS_DCN)
title('SS_corr_DCN', 'Interpreter', 'none')

ESN_Beautify_Plot(hFig, [8, 6])
end

%% function plot_SS_firing_histogram
function plot_burster_pauser
clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
SS_firing_rate_all = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.SS_firing_rate_all;

primSac_onset_train_data_logic_SS_all = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_SS_all;
% Definition-1: Saccade duration compared to baseline activity
% primSac_onset_activity = nanmean(primSac_onset_train_data_logic_SS_all(:,300:350), 2);
% primSac_onset_activity_change = (primSac_onset_activity*1000) - SS_firing_rate_all;
% Definition-2: Pre-saccade compare to post-saccade
primSac_onset_activity = nanmean(primSac_onset_train_data_logic_SS_all(:,300:400), 2);
primSac_onset_activity_prior = nanmean(primSac_onset_train_data_logic_SS_all(:,100:200), 2);
primSac_onset_activity_change = (primSac_onset_activity*1000) - (primSac_onset_activity_prior*1000);
% reshape and average
primSac_onset_activity_change = reshape(primSac_onset_activity_change, num_pCells, []);
primSac_onset_activity_change = nanmean(primSac_onset_activity_change, 2);

corrSac_onset_train_data_logic_SS_all = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.train_data_logic_SS_all;
% Definition-1: Saccade duration compared to baseline activity
% corrSac_onset_activity = nanmean(corrSac_onset_train_data_logic_SS_all(:,300:350), 2);
% corrSac_onset_activity_change = (corrSac_onset_activity*1000) - SS_firing_rate_all;
% Definition-2: Pre-saccade compare to post-saccade
corrSac_onset_activity = nanmean(corrSac_onset_train_data_logic_SS_all(:,300:400), 2);
corrSac_onset_activity_prior = nanmean(corrSac_onset_train_data_logic_SS_all(:,100:200), 2);
corrSac_onset_activity_change = (corrSac_onset_activity*1000) - (corrSac_onset_activity_prior*1000);
% reshape and average
corrSac_onset_activity_change = reshape(corrSac_onset_activity_change, num_pCells, []);
corrSac_onset_activity_change = nanmean(corrSac_onset_activity_change, 2);

activity_change = nanmean([primSac_onset_activity_change, corrSac_onset_activity_change], 2);

hFig = figure(103);
clf(hFig)
edges_ = -85:10:85;
subplot(1,3,1)
hold on
histogram(activity_change, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'm')
histogram(activity_change, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'm', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([0, 0], yl_, '-k', 'LineWidth', 1.5);
xlabel({'Change in activity', 'during saccade (Hz)'})
ylabel('Count (#)')
title('Prim+Corr')
subplot(1,3,2)
hold on
histogram(primSac_onset_activity_change, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'b')
histogram(primSac_onset_activity_change, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([0, 0], yl_, '-k', 'LineWidth', 1.5);
xlabel({'Change in activity', 'during saccade (Hz)'})
ylabel('Count (#)')
title('Prim')
subplot(1,3,3)
hold on
histogram(corrSac_onset_activity_change, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(corrSac_onset_activity_change, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([0, 0], yl_, '-k', 'LineWidth', 1.5);
xlabel({'Change in activity', 'during saccade (Hz)'})
ylabel('Count (#)')
title('Corr')

sgtitle(hFig, 'ALL Directions', 'Interpreter', 'none');
ESN_Beautify_Plot(hFig, [6, 4])

%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
SS_firing_rate_180 = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.SS_firing_rate_180;

primSac_onset_train_data_logic_SS_180 = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_SS_180;
% Definition-1: Saccade duration compared to baseline activity
% primSac_onset_activity = nanmean(primSac_onset_train_data_logic_SS_180(:,300:350), 2);
% primSac_onset_activity_change = (primSac_onset_activity*1000) - SS_firing_rate_180;
% Definition-2: Pre-saccade compare to post-saccade
primSac_onset_activity = nanmean(primSac_onset_train_data_logic_SS_180(:,300:400), 2);
primSac_onset_activity_prior = nanmean(primSac_onset_train_data_logic_SS_180(:,100:200), 2);
primSac_onset_activity_change = (primSac_onset_activity*1000) - (primSac_onset_activity_prior*1000);
% reshape and average
primSac_onset_activity_change = reshape(primSac_onset_activity_change, num_pCells, []);
primSac_onset_activity_change = nanmean(primSac_onset_activity_change, 2);

corrSac_onset_train_data_logic_SS_180 = RASTER_DATA_ALL_PCELL_TUNED.raster_data_corrSac_onset.train_data_logic_SS_180;
% Definition-1: Saccade duration compared to baseline activity
% corrSac_onset_activity = nanmean(corrSac_onset_train_data_logic_SS_180(:,300:350), 2);
% corrSac_onset_activity_change = (corrSac_onset_activity*1000) - SS_firing_rate_180;
% Definition-2: Pre-saccade compare to post-saccade
corrSac_onset_activity = nanmean(corrSac_onset_train_data_logic_SS_180(:,300:400), 2);
corrSac_onset_activity_prior = nanmean(corrSac_onset_train_data_logic_SS_180(:,100:200), 2);
corrSac_onset_activity_change = (corrSac_onset_activity*1000) - (corrSac_onset_activity_prior*1000);
% reshape and average
corrSac_onset_activity_change = reshape(corrSac_onset_activity_change, num_pCells, []);
corrSac_onset_activity_change = nanmean(corrSac_onset_activity_change, 2);

activity_change = nanmean([primSac_onset_activity_change, corrSac_onset_activity_change], 2);

hFig = figure(104);
clf(hFig)
edges_ = -95:10:95;
subplot(1,3,1)
hold on
histogram(activity_change, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'm')
histogram(activity_change, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'm', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([0, 0], yl_, '-k', 'LineWidth', 1.5);
xlabel({'Change in activity', 'during saccade (Hz)'})
ylabel('Count (#)')
title('Prim+Corr')
subplot(1,3,2)
hold on
histogram(primSac_onset_activity_change, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'b')
histogram(primSac_onset_activity_change, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([0, 0], yl_, '-k', 'LineWidth', 1.5);
xlabel({'Change in activity', 'during saccade (Hz)'})
ylabel('Count (#)')
title('Prim')
subplot(1,3,3)
hold on
histogram(corrSac_onset_activity_change, edges_,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(corrSac_onset_activity_change, edges_,  'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 2)
yl_ = ylim;
plot([0, 0], yl_, '-k', 'LineWidth', 1.5);
xlabel({'Change in activity', 'during saccade (Hz)'})
ylabel('Count (#)')
title('Corr')

sgtitle(hFig, 'CS+180 Direction', 'Interpreter', 'none');
ESN_Beautify_Plot(hFig, [6, 4])

%%
% extract burster/pauser cells and form a new ALL_PCELL_COMPRESSED_DATA
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);

primSac_onset_train_data_logic_SS_all = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_onset.train_data_logic_SS_all;
% Definition-2: Pre-saccade compare to post-saccade
primSac_onset_activity = nanmean(primSac_onset_train_data_logic_SS_all(:,250:400), 2);
primSac_onset_activity_prior = nanmean(primSac_onset_train_data_logic_SS_all(:,100:200), 2);
primSac_onset_activity_change = (primSac_onset_activity*1000) - (primSac_onset_activity_prior*1000);
% reshape and average
primSac_onset_activity_change = reshape(primSac_onset_activity_change, num_pCells, []);
primSac_onset_activity_change = nanmean(primSac_onset_activity_change, 2);

cell_idx_burster = find(primSac_onset_activity_change>0);
cell_idx_pauser = find(primSac_onset_activity_change<0);
cell_idx_interest = cell_idx_burster;

% ALL_PCELL_COMPRESSED_DATA_BKP = ALL_PCELL_COMPRESSED_DATA;
% ALL_PCELL_COMPRESSED_DATA = ALL_PCELL_COMPRESSED_DATA_BKP;
clearvars ALL_PCELL_COMPRESSED_DATA;
for counter_cell = 1 : length(cell_idx_interest)
    ALL_PCELL_COMPRESSED_DATA(counter_cell) = ALL_PCELL_COMPRESSED_DATA_BKP(cell_idx_interest(counter_cell));
end

end

%% function plot_CS_on_two_bump
function plot_CS_on_two_bump
clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
%%
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
cue_present_train_data_logic_CS_000 = RASTER_DATA_ALL_PCELL_TUNED.raster_data_cue_present.train_data_logic_CS_000;
cue_present_smooth = smoothdata(cue_present_train_data_logic_CS_000, 2, 'movmean', 21, 'omitnan');
% cue_present_smooth = cue_present_train_data_logic_CS_000;
hFig = figure(107);
clf(hFig)
subplot(2,2,1)
hold on
time_cue_present = (-50:1:250);
range_ = (250:1:550);
cue_present_smooth = (cue_present_smooth(:,range_).*1000);
plot_mean = nanmean(cue_present_smooth);
plot_stdv = nanstd(cue_present_smooth)./sqrt(num_pCells);
plot_stdv_p = plot_mean + plot_stdv;
plot_stdv_m = plot_mean - plot_stdv;
plot(time_cue_present, plot_stdv_m, 'k')
plot(time_cue_present, plot_stdv_p, 'k')
plot(time_cue_present, plot_mean, 'k')
yl_ = ylim;
plot([45 45], yl_, 'k');
plot([95 95], yl_, 'k');
plot([145 145], yl_, 'k');
xlabel('Time for cue (ms)')
ylabel('CS firing (Hz)')

subplot(2,2,2)
hold on
peak_1 = mean(cue_present_train_data_logic_CS_000(:,(351:400)-5), 2)*1000;
peak_2 = mean(cue_present_train_data_logic_CS_000(:,(401:450)-5), 2)*1000;

x = peak_1(:);
y = peak_2(:);
slope = polyfit(x,y,1);
yfit = polyval(slope,sort(x));
mdl = fitlm(x,y);
p_val = mdl.Coefficients.pValue;
slope = mdl.Coefficients.Estimate;

plot(peak_1, peak_2, '.k');
plot(sort(x), yfit, '-k');
xlabel('Avg. firing 45-95ms')
ylabel('Avg. firing 95-145ms')
title({['slope:' mat2str(slope,3)], ['p-value:' mat2str(p_val,3)]})

clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
num_pCells = length(ALL_PCELL_COMPRESSED_DATA);
primSac_offset_train_data_logic_CS_000 = RASTER_DATA_ALL_PCELL_TUNED.raster_data_primSac_offset.train_data_logic_CS_000;
primSac_offset_smooth = smoothdata(primSac_offset_train_data_logic_CS_000, 2, 'movmean', 21, 'omitnan');
% cue_present_smooth = cue_present_train_data_logic_CS_000;
hFig = figure(107);
subplot(2,2,3)
hold on
time_primSac_offset = (-50:1:250);
range_ = (250:1:550);
primSac_offset_smooth = (primSac_offset_smooth(:,range_).*1000);
plot_mean = nanmean(primSac_offset_smooth);
plot_stdv = nanstd(primSac_offset_smooth)./sqrt(num_pCells);
plot_stdv_p = plot_mean + plot_stdv;
plot_stdv_m = plot_mean - plot_stdv;
plot(time_primSac_offset, plot_stdv_m, 'k')
plot(time_primSac_offset, plot_stdv_p, 'k')
plot(time_primSac_offset, plot_mean, 'k')
yl_ = ylim;
plot([30 30], yl_, 'k');
plot([80 80], yl_, 'k');
plot([130 130], yl_, 'k');
xlabel('Time for primSac offset (ms)')
ylabel('CS firing (Hz)')

subplot(2,2,4)
hold on
peak_1 = mean(primSac_offset_train_data_logic_CS_000(:,(351:400)-20), 2)*1000;
peak_2 = mean(primSac_offset_train_data_logic_CS_000(:,(401:450)-20), 2)*1000;

x = peak_1(:);
y = peak_2(:);
slope = polyfit(x,y,1);
yfit = polyval(slope,sort(x));
mdl = fitlm(x,y);
p_val = mdl.Coefficients.pValue;
slope = mdl.Coefficients.Estimate;

plot(peak_1, peak_2, '.k');
plot(sort(x), yfit, '-k');
xlabel('Avg. firing 30-80ms')
ylabel('Avg. firing 80-130ms')
title({['slope:' mat2str(slope,3)], ['p-value:' mat2str(p_val,3)]})

ESN_Beautify_Plot(hFig, [6, 4])
end

%% function plot_bootstrap_population_CS000_CS180
function plot_bootstrap_population_CS000_CS180
clearvars -except path_data_monkey_sorted pCell_list ALL_PCELL_COMPRESSED_DATA RASTER_DATA_ALL_PCELL_TUNED
%%
field_name = 'raster_data_primSac_onset';
train_data_dir = '_180';
numTrial = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['numTrial' train_data_dir]);
SS_firing_rate = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['SS_firing_rate' train_data_dir]);
train_data_logic_SS = RASTER_DATA_ALL_PCELL_TUNED.(field_name).(['train_data_logic_SS' train_data_dir]);
train_data_logic_SS = train_data_logic_SS - (SS_firing_rate ./ 1000);

nBoot = 2000; %number of bootstraps
num_time_points = size(train_data_logic_SS, 2);
variable_data_mean = nan(1, num_time_points);
variable_data_stdv_p = nan(1, num_time_points);
variable_data_stdv_m = nan(1, num_time_points);
tic
for counter_time = 1 : num_time_points
    data_ = train_data_logic_SS(:,counter_time);
    [bci,bmeans] = bootci(nBoot,{@mean,data_},'alpha',.05,'type','per'); %90 confidence interval
    bmu = mean(bmeans); % Compute bootstrap sample mean
    variable_data_mean(1, counter_time) = bmu;
    variable_data_stdv_p(1, counter_time) = bci(2);
    variable_data_stdv_m(1, counter_time) = bci(1);
end
toc

hFig = figure(108);
clf(hFig)
hold on
plot(ESN_smooth( variable_data_stdv_m ).*1000, 'k')
plot(ESN_smooth( variable_data_stdv_p ).*1000, 'k')
plot(ESN_smooth( variable_data_mean ).*1000, 'k')

%%
num_perm = 5000;
num_samples = 50; % num_pCells
length_span = size(train_data_logic_SS, 2);
train_data_logic_SS_perm = nan(num_perm, length_span);
for counter_perm = 1 : num_perm
    inds_perm = randi(num_pCells,[num_samples 1]);
    train_data_logic_SS_ = train_data_logic_SS(inds_perm, :);
    train_data_logic_SS_(isnan(train_data_logic_SS_)) = 0;
    numTrial_ = numTrial(inds_perm, :);
    firing_SS_mean = numTrial_' * train_data_logic_SS_ ./ nansum(numTrial_) .* 1000;
    train_data_logic_SS_perm(counter_perm, :) = firing_SS_mean;
end
variable_data_mean = ESN_smooth( nanmean(train_data_logic_SS_perm, 1) );
variable_data_stdv = ESN_smooth( nanstd(train_data_logic_SS_perm, 0, 1) ) ./ sqrt(num_samples);

end

%% function plot_CS_bandpower
function plot_cs_bandpower()
%%
psort_data = Psort_read_psort;
%% cs_waveform
ch_data = psort_data.topLevel_data.ch_data;
sample_rate = psort_data.topLevel_data.sample_rate;
win_len_before = 0.5e-3; % 30e-3; % 
win_len_after = 10e-3; % 30e-3; % 
% cs_waveform
spike_bool = psort_data.topLevel_data.cs_index;
[cs_waveform, ~] = extract_waveform(ch_data, spike_bool, sample_rate, win_len_before, win_len_after);
cs_waveform = cs_waveform - median(cs_waveform);
% ss_waveform
spike_bool = psort_data.topLevel_data.ss_index;
[ss_waveform, ~] = extract_waveform(ch_data, spike_bool, sample_rate, win_len_before, win_len_after);
% ss_waveform = ss_waveform - median(ss_waveform);
%% cs_power
pband = bandpower(cs_waveform',sample_rate,[0 100]);
ptot  = bandpower(cs_waveform',sample_rate,[0 5000]);
cs_power = 100*(pband./ptot);
pband = bandpower(ss_waveform',sample_rate,[0 100]);
ptot  = bandpower(ss_waveform',sample_rate,[0 5000]);
ss_power = 100*(pband./ptot);
%% plot power
hFig = figure(1);
clf(hFig);
subplot(3,2,1)
plot(cs_power, '.k');
subplot(3,2,2)
plot(ss_power, '.k');
subplot(3,2,3)
plot(detrend(cs_power), '.k');
subplot(3,2,4)
plot(detrend(ss_power), '.k');
subplot(3,2,5)
hold on
histogram((cs_power),  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram((cs_power),  'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 2)
subplot(3,2,6)
hold on
histogram((ss_power),  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'b')
histogram((ss_power),  'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'FaceColor', 'none', 'linewidth', 2)

ESN_Beautify_Plot(hFig, [8,6])

%% plot spect
win_len_before_ = 30e-3; % 0.5e-3; % 
win_len_after_ = 30e-3; % 10e-3; % 
% cs_waveform
spike_bool_ = psort_data.topLevel_data.cs_index;
[cs_waveform_, cs_span_] = extract_waveform(ch_data, spike_bool_, sample_rate, win_len_before_, win_len_after_);
% cs_waveform_ = cs_waveform_ - median(cs_waveform_);
%%
idx = 1502;
signal = cs_waveform_(idx,:)';
[s,f,t] = spectrogram(signal,32,30,2^12,sample_rate);
t = (t - mean(t) + mean(mean(cs_span_))) * 1000;

hFig = figure(2);
clf(hFig);
subplot(4,1,1)
plot(mean(cs_span_)*1000, signal, '-k')
ylim([-2000 2500])
xlim([-30 30])
ylabel('Waveform (uV)')

subplot(4,1,[2 3])
surf(t, f, abs(s));
shading interp
view([0 0 90])
ylim([0 5000])
xlim([-30 30])
ylabel('Frequency (Hz)')

subplot(4,1,4)
surf(t, f, abs(s));
shading interp
view([0 0 90])
ylim([0 100])
xlim([-30 30])
xlabel('Time (ms)')

end

%% function extract_waveform
function [waveform, span] = extract_waveform(data, spike_bool, sample_rate, win_len_before, win_len_after)
spike_bool = logical(spike_bool); spike_bool(1) = false; spike_bool(end) = false;
spike_int = find(spike_bool);
win_len_before_int = round(double(win_len_before) * double(sample_rate));
win_len_after_int  = round(double(win_len_after)  * double(sample_rate));
span_int = (-win_len_before_int : 1 : win_len_after_int);
num_row = length(spike_int);
num_col = length(span_int);
spike_int = repmat(spike_int(:), 1, num_col);
span_int  = repmat(span_int(:)', num_row, 1);

ind = spike_int + span_int;
ind(ind<1) = 1;
ind(ind>length(data)) = length(data);
waveform = data(ind);
span = span_int / double(sample_rate);
end

%% copy plot_data file
function copy_plot_data()
%% build and save ALL_PCELL_COMPRESSED_DATA
path_data_monkey_sorted = uigetdir;

% build pCell_list, this is a hard coded cell with the id of all of the pCells and the bundles
pCell_list_1 = build_pCell_list_Mirza_pre201906();
pCell_list_2 = build_pCell_list_Mirza_post201906();
pCell_list_3 = build_pCell_list_Ramon();
pCell_list_4 = build_pCell_list_Mirza_post202011();
pCell_list = vertcat(pCell_list_1, pCell_list_2, pCell_list_3, pCell_list_4);
%% init vars
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ... ']);
    num_recording = sum(pCell_list_isstr(counter_pCell, :));
    plot_data_cell = struct();
    counter_recording = 1;
    %% build plot_data address
    file_name = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN';
    year_ = file_name(1:2);
    month_ = file_name(3:4);
    day_ = file_name(5:6);
    hour_ = file_name(8:9);
    minute_ = file_name(10:11);
    second_ = file_name(12:13);
    subFolder_month = ['20' year_ '-' month_ filesep];
    subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
    subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
    if num_recording < 2
        subFolder_figs = ['analyzed_figs' filesep];
        file_name = [file_name '_plot_data' '.mat']; % '190423_142023_01_sorted_ESN_plot_data';
    else
        subFolder_figs = ['bundle_figs' filesep];
        file_name = [file_name '_combine_' num2str(num_recording) '_plot_data' '.mat']; % '190326_180641_03_sorted_RSh_combine_2_plot_data'
    end
    file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
    file_full_path = [file_path file_name];
    file_full_destination = [path_data_monkey_sorted 'pCell_properties' filesep file_name];
    
    copyfile(file_full_path, file_full_destination)
    fprintf(' --> Completed. \n')
end
end