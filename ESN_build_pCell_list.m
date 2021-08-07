%% function ESN_build_pCell_list
function pCell_list = ESN_build_pCell_list(flag_pair_list)
%% Handle inputs
if nargin < 1
    flag_pair_list = false;
end

%% build pCell_list, this is a hard coded cell with the id of all of the pCells and the bundles
if ~flag_pair_list
%% build
pCell_list_1 = build_pCell_list_Mirza_pre201906();
pCell_list_2 = build_pCell_list_Mirza_post201906();
pCell_list_3 = build_pCell_list_Mirza_post202011();
pCell_list_4 = build_pCell_list_Ramon();
pCell_list = vertcat(pCell_list_1, pCell_list_2, pCell_list_3, pCell_list_4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build pair list
if flag_pair_list
%% build
pair_list_full = build_pair_list_full();
pCell_list = pair_list_full;
end
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

pCell_list_Mirza_post201906(11, 1:1) = {'190819_154221_01_sorted_MAF'};
pCell_list_Mirza_post201906(12, 1:1) = {'190819_154221_03_sorted_MAF'};
pCell_list_Mirza_post201906(13, 1:1) = {'190820_164746_01_sorted_MAF'};
pCell_list_Mirza_post201906(14, 1:3) = {'190823_160211_01_sorted_MAF', '190823_162607_01_sorted_MAF', '190823_165201_01_sorted_MAF'};
pCell_list_Mirza_post201906(15, 1:3) = {'190823_160211_02_sorted_MAF', '190823_162607_02_sorted_MAF', '190823_165201_02_sorted_MAF'};
pCell_list_Mirza_post201906(16, 1:1) = {'190823_142416_03_sorted_MAF'};
pCell_list_Mirza_post201906(17, 1:2) = {'190918_115126_03_sorted_RSh', '190918_122153_03_sorted_RSh'};
pCell_list_Mirza_post201906(18, 1:2) = {'190918_115126_04_sorted_RSh', '190918_122153_04_sorted_RSh'};
pCell_list_Mirza_post201906(19, 1:2) = {'190919_095155_02_sorted_RSh', '190919_100229_02_sorted_RSh'};
pCell_list_Mirza_post201906(20, 1:3) = {'190923_143725_02_sorted_PGH', '190923_150658_02_sorted_PGH', '190923_154446_02_sorted_PGH'};

pCell_list_Mirza_post201906(21, 1:3) = {'190925_145115_02_sorted_PGH', '190925_152228_02_sorted_PGH', '190925_155449_02_sorted_PGH'};
pCell_list_Mirza_post201906(22, 1:1) = {'190927_135905_01_sorted_RSh'};
pCell_list_Mirza_post201906(23, 1:1) = {'190927_135905_03_sorted_RSh'};


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

pCell_list_Ramon(11, 1:1) = {'190919_151958_03_sorted_RSh'};
pCell_list_Ramon(12, 1:1) = {'190919_151958_04_sorted_RSh'};
pCell_list_Ramon(13, 1:1) = {'190919_154206_02_sorted_ESN'};
pCell_list_Ramon(14, 1:2) = {'190925_113310_02_sorted_ESN', '190925_121201_02_sorted_ESN'};
pCell_list_Ramon(15, 1:1) = {'191010_122539_05_sorted_ESN'};
pCell_list_Ramon(16, 1:1) = {'191010_125057_02_sorted_ESN'};
pCell_list_Ramon(17, 1:3) = {'191101_130349_03_sorted_ESN', '191101_133125_03_sorted_ESN', '191101_140646_03_sorted_ESN'};
pCell_list_Ramon(18, 1:3) = {'191101_130349_04_sorted_ESN', '191101_133125_04_sorted_ESN', '191101_140646_04_sorted_ESN'};
pCell_list_Ramon(19, 1:3) = {'191101_130349_05_sorted_ESN', '191101_133125_05_sorted_ESN', '191101_140646_05_sorted_ESN'};
pCell_list_Ramon(20, 1:2) = {'191104_123450_02_sorted_ESN', '191104_125346_02_sorted_ESN'};

pCell_list_Ramon(21, 1:2) = {'191104_123450_07_sorted_RSh', '191104_125346_07_sorted_RSh'};
pCell_list_Ramon(22, 1:1) = {'191104_131449_03_sorted_ESN'};
pCell_list_Ramon(23, 1:2) = {'191108_113612_04_sorted_ESN', '191108_114524_04_sorted_ESN'};
pCell_list_Ramon(24, 1:1) = {'191108_124201_03_sorted_ESN'};
pCell_list_Ramon(25, 1:1) = {'191108_124201_04_sorted_ESN'};
pCell_list_Ramon(26, 1:1) = {'191108_135053_07_sorted_ESN'};
pCell_list_Ramon(27, 1:4) = {'191111_131030_04_sorted_ESN', '191111_132846_04_sorted_ESN', '191111_135712_04_sorted_ESN', '191111_142741_04_sorted_ESN'};
pCell_list_Ramon(28, 1:1) = {'191111_131030_04_2_sorted_ESN'};
pCell_list_Ramon(29, 1:2) = {'191111_135712_03_sorted_RSh','191111_142741_03_sorted_RSh'};
pCell_list_Ramon(30, 1:8) = {'191115_115149_04_sorted_RSh', '191115_120043_04_sorted_RSh', '191115_120257_04_sorted_RSh', '191115_122302_04_sorted_RSh', '191115_122724_04_sorted_RSh', '191115_124011_04_sorted_RSh', '191115_125241_04_sorted_RSh', '191115_125423_04_sorted_RSh'};

pCell_list_Ramon(31, 1:3) = {'191118_143441_05_sorted_RSh', '191118_145837_05_sorted_RSh', '191118_154644_05_sorted_RSh'};
pCell_list_Ramon(32, 1:3) = {'191118_145837_01_sorted_RSh', '191118_152411_01_sorted_RSh', '191118_154644_01_sorted_RSh'};
pCell_list_Ramon(33, 1:3) = {'191120_141613_02_sorted_RSh', '191120_144254_02_sorted_RSh', '191120_145548_02_sorted_RSh'};
pCell_list_Ramon(34, 1:1) = {'191204_115824_03_sorted_RSh'};
pCell_list_Ramon(35, 1:1) = {'191204_131009_01_sorted_RSh'};
pCell_list_Ramon(36, 1:1) = {'191209_141244_03_sorted_JSP'};
pCell_list_Ramon(37, 1:1) = {'191209_141244_04_sorted_RSh'};
pCell_list_Ramon(38, 1:2) = {'200109_110904_04_sorted_RSh', '200109_112601_04_sorted_RSh'};
pCell_list_Ramon(39, 1:1) = {'200109_125642_03_sorted_RSh'};
pCell_list_Ramon(40, 1:3) = {'200122_125828_07_sorted_RSh', '200122_133020_07_sorted_RSh', '200122_135404_07_sorted_RSh'};

pCell_list_Ramon(41, 1:2) = {'200122_133020_04_sorted_RSh', '200122_135404_04_sorted_RSh'};
pCell_list_Ramon(42, 1:4) = {'200123_124808_07_sorted_PGH', '200123_131422_07_sorted_PGH', '200123_133835_07_sorted_PGH', '200123_140924_07_sorted_PGH'};
pCell_list_Ramon(43, 1:4) = {'200128_151251_04_sorted_RSh', '200128_151937_06_sorted_RSh', '200128_153045_06_sorted_RSh', '200128_154444_01_sorted_RSh'};
pCell_list_Ramon(44, 1:1) = {'200128_163856_04_sorted_RSh'};
pCell_list_Ramon(45, 1:4) = {'200129_094755_04_sorted_PGH', '200129_124934_04_sorted_PGH', '200129_131018_04_sorted_PGH', '200129_133230_04_sorted_PGH'};
pCell_list_Ramon(46, 1:3) = {'200204_162413_04_sorted_RSh', '200204_163427_04_sorted_RSh', '200204_164913_03_sorted_RSh'};
pCell_list_Ramon(47, 1:1) = {'200211_162022_03_sorted_RSh'};
pCell_list_Ramon(48, 1:5) = {'200303_153615_03_sorted_RSh', '200303_155436_03_sorted_RSh', '200303_163555_03_sorted_RSh', '200303_165850_03_sorted_RSh', '200303_172913_03_sorted_RSh'};
pCell_list_Ramon(49, 1:3) = {'200814_113452_05_sorted_RSh', '200814_114000_05_sorted_RSh', '200814_114425_05_sorted_RSh'};
pCell_list_Ramon(50, 1:3) = {'200914_122203_04_sorted_RSh', '200914_124954_04_sorted_RSh', '200914_132356_04_sorted_RSh'};

pCell_list_Ramon(51, 1:3) = {'200916_130045_02_sorted_PGH', '200916_133908_02_sorted_PGH', '200916_142429_02_sorted_PGH'};
pCell_list_Ramon(52, 1:1) = {'200923_133039_03_sorted_PGH'};
pCell_list_Ramon(53, 1:4) = {'201007_132359_04_sorted_PGH', '201007_133319_02_sorted_PGH', '201007_134353_02_sorted_PGH', '201007_135402_04_sorted_PGH'};

end

%% function build_pair_list_full
function pair_list_full = build_pair_list_full()

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_list_full_Mirza_pre201906 = cell(0,10);
pair_list_full_Mirza_pre201906( 1, 1:2) = {'190326_180641_03_sorted_RSh', '190326_182755_03_sorted_RSh'};
pair_list_full_Mirza_pre201906( 2, 1:2) = {'190326_180641_04_sorted_RSh', '190326_182755_04_sorted_RSh'};
pair_list_full_Mirza_pre201906( 3, 1:3) = {'190327_135803_05_sorted_RSh', '190327_143005_05_sorted_RSh', '190327_150203_05_sorted_RSh'};
pair_list_full_Mirza_pre201906( 4, 1:3) = {'190327_135803_06_sorted_RSh', '190327_143005_06_sorted_RSh', '190327_150203_06_sorted_RSh'};
pair_list_full_Mirza_pre201906( 5, 1:1) = {'190329_133538_05_sorted_RSh'};
pair_list_full_Mirza_pre201906( 6, 1:1) = {'190329_133538_03_sorted_RSh'};
pair_list_full_Mirza_pre201906( 7, 1:1) = {'190329_150427_03_sorted_RSh'};
pair_list_full_Mirza_pre201906( 8, 1:1) = {'190329_150427_07_sorted_RSh'};
pair_list_full_Mirza_pre201906( 9, 1:1) = {'190409_150223_03_sorted_ESN'};
pair_list_full_Mirza_pre201906(10, 1:1) = {'190409_150223_07_sorted_ESN'};

pair_list_full_Mirza_pre201906(11, 1:1) = {'190409_152524_03_sorted_ESN'};
pair_list_full_Mirza_pre201906(12, 1:1) = {'190409_152524_07_sorted_ESN'};
pair_list_full_Mirza_pre201906(13, 1:1) = {'190411_130256_03_sorted_RSh'};
pair_list_full_Mirza_pre201906(14, 1:1) = {'190411_130256_07_sorted_ESN'};
pair_list_full_Mirza_pre201906(15, 1:5) = {'190412_115230_05_sorted_RSh', '190412_120702_05_sorted_RSh', '190412_122025_05_sorted_RSh', '190412_123356_05_sorted_RSh', '190412_125020_05_sorted_RSh'};
pair_list_full_Mirza_pre201906(16, 1:5) = {'190412_115230_01_sorted_RSh', '190412_120702_03_sorted_RSh', '190412_122025_03_sorted_RSh', '190412_123356_03_sorted_RSh', '190412_125020_01_sorted_RSh'};
pair_list_full_Mirza_pre201906(17, 1:1) = {'190422_144948_05_sorted_RSh'};
pair_list_full_Mirza_pre201906(18, 1:1) = {'190422_144948_07_sorted_RSh'};
pair_list_full_Mirza_pre201906(19, 1:2) = {'190426_120355_05_sorted_RSh', '190426_122359_05_sorted_RSh'};
pair_list_full_Mirza_pre201906(20, 1:2) = {'190426_120355_06_sorted_RSh', '190426_122359_06_sorted_RSh'};

pair_list_full_Mirza_pre201906(21, 1:3) = {'190426_125234_07_sorted_RSh', '190426_131912_07_sorted_RSh', '190426_134533_07_sorted_RSh'};
pair_list_full_Mirza_pre201906(22, 1:3) = {'190426_125234_01_sorted_RSh', '190426_131912_01_sorted_RSh', '190426_134533_01_sorted_RSh'};
pair_list_full_Mirza_pre201906(23, 1:2) = {'190429_135343_05_sorted_RSh', '190429_141529_05_sorted_RSh'};
pair_list_full_Mirza_pre201906(24, 1:2) = {'190429_135343_06_sorted_RSh', '190429_141529_06_sorted_RSh'};
pair_list_full_Mirza_pre201906(25, 1:3) = {'190429_142903_02_sorted_RSh', '190429_145043_02_sorted_JSP', '190429_150541_02_sorted_JSP'};
pair_list_full_Mirza_pre201906(26, 1:3) = {'190429_142903_03_sorted_RSh', '190429_145043_03_sorted_RSh', '190429_150541_03_sorted_RSh'};
pair_list_full_Mirza_pre201906(27, 1:3) = {'190429_142903_02_sorted_RSh', '190429_145043_02_sorted_JSP', '190429_150541_02_sorted_JSP'};
pair_list_full_Mirza_pre201906(28, 1:3) = {'190429_142903_07_sorted_RSh', '190429_145043_07_sorted_RSh', '190429_150541_07_sorted_RSh'};
pair_list_full_Mirza_pre201906(29, 1:3) = {'190429_142903_03_sorted_RSh', '190429_145043_03_sorted_RSh', '190429_150541_03_sorted_RSh'};
pair_list_full_Mirza_pre201906(30, 1:3) = {'190429_142903_07_sorted_RSh', '190429_145043_07_sorted_RSh', '190429_150541_07_sorted_RSh'};

pair_list_full_Mirza_pre201906(31, 1:2) = {'190516_144904_02_sorted_RSh', '190516_150909_02_sorted_RSh'};
pair_list_full_Mirza_pre201906(32, 1:2) = {'190516_144904_07_sorted_RSh', '190516_150909_07_sorted_RSh'};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_list_full_Mirza_post201906 = cell(0,10);
pair_list_full_Mirza_post201906( 1, 1:1) = {'190722_154305_03_sorted_RSh'};
pair_list_full_Mirza_post201906( 2, 1:1) = {'190722_154305_07_sorted_PGH'};
pair_list_full_Mirza_post201906( 3, 1:1) = {'190819_154221_01_sorted_MAF'};
pair_list_full_Mirza_post201906( 4, 1:1) = {'190819_154221_03_sorted_MAF'};
pair_list_full_Mirza_post201906( 5, 1:3) = {'190823_160211_01_sorted_MAF', '190823_162607_01_sorted_MAF', '190823_165201_01_sorted_MAF'};
pair_list_full_Mirza_post201906( 6, 1:3) = {'190823_160211_02_sorted_MAF', '190823_162607_02_sorted_MAF', '190823_165201_02_sorted_MAF'};
pair_list_full_Mirza_post201906( 7, 1:2) = {'190918_115126_03_sorted_RSh', '190918_122153_03_sorted_RSh'};
pair_list_full_Mirza_post201906( 8, 1:2) = {'190918_115126_04_sorted_RSh', '190918_122153_04_sorted_RSh'};
pair_list_full_Mirza_post201906( 9, 1:1) = {'190927_135905_01_sorted_RSh'};
pair_list_full_Mirza_post201906(10, 1:1) = {'190927_135905_03_sorted_RSh'};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_list_full_Mirza_post202011 = cell(0,10);
pair_list_full_Mirza_post202011( 1, 1:1) = {'201223_131532_07_sorted_RSh'};
pair_list_full_Mirza_post202011( 2, 1:1) = {'201223_131532_03_sorted_RSh'};
pair_list_full_Mirza_post202011( 3, 1:2) = {'210104_151914_01_sorted_RSh', '210104_153305_01_sorted_RSh'};
pair_list_full_Mirza_post202011( 4, 1:2) = {'210104_151914_02_sorted_RSh', '210104_153305_02_sorted_RSh'};
pair_list_full_Mirza_post202011( 5, 1:5) = {'210106_141642_03_sorted_RSh', '210106_142808_03_sorted_RSh', '210106_145826_03_sorted_RSh', '210106_150439_03_sorted_RSh', '210106_153603_03_sorted_RSh'};
pair_list_full_Mirza_post202011( 6, 1:5) = {'210106_141642_05_sorted_RSh', '210106_142808_05_sorted_RSh', '210106_145826_05_sorted_RSh', '210106_150439_05_sorted_RSh', '210106_153603_05_sorted_RSh'};
pair_list_full_Mirza_post202011( 7, 1:1) = {'210107_154500_03_sorted_RSh'};
pair_list_full_Mirza_post202011( 8, 1:1) = {'210107_154500_06_sorted_RSh'};
pair_list_full_Mirza_post202011( 9, 1:1) = {'210107_154500_03_sorted_RSh'};
pair_list_full_Mirza_post202011(10, 1:1) = {'210107_154500_07_sorted_RSh'};

pair_list_full_Mirza_post202011(11, 1:1) = {'210107_154500_06_sorted_RSh'};
pair_list_full_Mirza_post202011(12, 1:1) = {'210107_154500_07_sorted_RSh'};
pair_list_full_Mirza_post202011(13, 1:2) = {'210108_151142_06_sorted_RSh', '210108_153618_06_sorted_RSh'};
pair_list_full_Mirza_post202011(14, 1:2) = {'210108_151142_07_sorted_RSh', '210108_153618_07_sorted_RSh'};
pair_list_full_Mirza_post202011(15, 1:2) = {'210122_144322_05_sorted_RSh', '210122_145740_05_sorted_RSh'};
pair_list_full_Mirza_post202011(16, 1:2) = {'210122_144322_06_sorted_RSh', '210122_145740_06_sorted_RSh'};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pair_list_full_Ramon = cell(0,10);
pair_list_full_Ramon( 1, 1:2) = {'190830_114101_02_sorted_ESN', '190830_114850_02_sorted_ESN'};
pair_list_full_Ramon( 2, 1:2) = {'190830_114101_04_sorted_ESN', '190830_114850_04_sorted_ESN'};
pair_list_full_Ramon( 3, 1:2) = {'190903_160127_04_sorted_ESN', '190903_163920_04_sorted_ESN'};
pair_list_full_Ramon( 4, 1:2) = {'190903_160127_01_sorted_ESN', '190903_163920_01_sorted_ESN'};
pair_list_full_Ramon( 5, 1:1) = {'190919_151958_03_sorted_RSh'};
pair_list_full_Ramon( 6, 1:1) = {'190919_151958_04_sorted_RSh'};
pair_list_full_Ramon( 7, 1:3) = {'191101_130349_03_sorted_ESN', '191101_133125_03_sorted_ESN', '191101_140646_03_sorted_ESN'};
pair_list_full_Ramon( 8, 1:3) = {'191101_130349_05_sorted_ESN', '191101_133125_05_sorted_ESN', '191101_140646_05_sorted_ESN'};
pair_list_full_Ramon( 9, 1:2) = {'191101_133125_03_sorted_ESN', '191101_140646_03_sorted_ESN'};
pair_list_full_Ramon(10, 1:2) = {'191101_133125_04_sorted_ESN', '191101_140646_04_sorted_ESN'};

pair_list_full_Ramon(11, 1:2) = {'191101_133125_04_sorted_ESN', '191101_140646_04_sorted_ESN'};
pair_list_full_Ramon(12, 1:2) = {'191101_133125_05_sorted_ESN', '191101_140646_05_sorted_ESN'};
pair_list_full_Ramon(13, 1:2) = {'191104_123450_02_sorted_ESN', '191104_125346_02_sorted_ESN'};
pair_list_full_Ramon(14, 1:2) = {'191104_123450_07_sorted_RSh', '191104_125346_07_sorted_RSh'};
pair_list_full_Ramon(15, 1:1) = {'191108_124201_03_sorted_ESN'};
pair_list_full_Ramon(16, 1:1) = {'191108_124201_04_sorted_ESN'};
pair_list_full_Ramon(17, 1:1) = {'191111_131030_04_sorted_ESN'};
pair_list_full_Ramon(18, 1:1) = {'191111_131030_04_2_sorted_ESN'};
pair_list_full_Ramon(19, 1:2) = {'191111_135712_03_sorted_RSh', '191111_142741_03_sorted_RSh'};
pair_list_full_Ramon(20, 1:2) = {'191111_135712_04_sorted_ESN', '191111_142741_04_sorted_ESN'};

pair_list_full_Ramon(21, 1:2) = {'191118_145837_05_sorted_RSh', '191118_154644_05_sorted_RSh'};
pair_list_full_Ramon(22, 1:2) = {'191118_145837_01_sorted_RSh', '191118_154644_01_sorted_RSh'};
pair_list_full_Ramon(23, 1:1) = {'191209_141244_03_sorted_JSP'};
pair_list_full_Ramon(24, 1:1) = {'191209_141244_04_sorted_RSh'};
pair_list_full_Ramon(25, 1:1) = {'200122_133020_07_sorted_RSh'};
pair_list_full_Ramon(26, 1:1) = {'200122_133020_04_sorted_RSh'};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pair_list_full_Jay = cell(0,10);
% pair_list_full_Jay( 1, 1:1) = {'191126_124046_02_sorted_JSP'};
% pair_list_full_Jay( 2, 1:1) = {'191126_124046_05_sorted_RSh'};
% pair_list_full_Jay( 3, 1:2) = {'200102_170006_03_sorted_RSh','200102_172224_03_sorted_RSh'};
% pair_list_full_Jay( 4, 1:2) = {'200102_170006_04_sorted_RSh','200102_172224_04_sorted_RSh'};
% pair_list_full_Jay( 5, 1:1) = {'200902_141429_01_sorted_RSh'};
% pair_list_full_Jay( 6, 1:1) = {'200902_141429_03_sorted_RSh'};
% pair_list_full_Jay( 7, 1:4) = {'200930_170108_01_sorted_JSP','200930_172823_01_sorted_JSP','200930_180153_01_sorted_JSP','200930_181219_01_sorted_JSP'};
% pair_list_full_Jay( 8, 1:4) = {'200930_170108_06_sorted_RSh','200930_172823_06_sorted_RSh','200930_180153_06_sorted_RSh','200930_181219_06_sorted_RSh'};
% pair_list_full_Jay( 9, 1:2) = {'210303_140235_05_sorted_RSh','210303_142154_05_sorted_RSh'};
% pair_list_full_Jay(10, 1:2) = {'210303_140235_07_sorted_RSh','210303_142154_07_sorted_RSh'};

% pair_list_full_Jay(11, 1:3) = {'210322_125218_04_sorted_RSh','210322_132133_04_sorted_RSh','210322_135448_04_sorted_RSh'};
% pair_list_full_Jay(12, 1:3) = {'210322_125218_06_sorted_RSh','210322_132133_06_sorted_RSh','210322_135448_06_sorted_RSh'};
% pair_list_full_Jay(13, 1:2) = {'210322_132133_03_sorted_RSh','210322_135448_03_sorted_RSh'};
% pair_list_full_Jay(14, 1:2) = {'210322_132133_04_sorted_RSh','210322_135448_04_sorted_RSh'};

pair_list_full = vertcat(pair_list_full_Mirza_pre201906, pair_list_full_Mirza_post201906, pair_list_full_Mirza_post202011, pair_list_full_Ramon);
% pair_list_full = vertcat(pair_list_full_Mirza_pre201906, pair_list_full_Mirza_post201906, pair_list_full_Mirza_post202011, pair_list_full_Ramon, pair_list_full_Jay);
end

