function ESN_Plot_Bkg(Y_Min_Max, Section_Values, Patch_Values, Colors, Block_trials)
hold on
if nargin < 5
    Block_trials = [];
end
Y_Min_Max = Y_Min_Max(:);
Y_Min = Y_Min_Max(1);
Y_Max = Y_Min_Max(2);

Section_Values = Section_Values(:);
Patch_Values = Patch_Values(:);
if(~isempty(Block_trials))
    for counter = 1 : length(Block_trials)
        patch([Block_trials(counter);Block_trials(counter);NaN],[Y_Min;Y_Max;NaN],'w','linewidth',1,'EdgeAlpha',0.2,'EdgeColor','k', 'lineStyle', '--');
    end
end
if(~isempty(Patch_Values))
    for counter = 1 : length(Patch_Values)
        patch([Section_Values(1);Section_Values(end);NaN],[Patch_Values(counter);Patch_Values(counter);NaN],'w','linewidth',2,'EdgeAlpha',0.6,'EdgeColor','k', 'lineStyle', '--');
    end
end

for counter = 1 : length(Section_Values)-1
    fill([Section_Values(counter) Section_Values(counter+1) Section_Values(counter+1) Section_Values(counter) Section_Values(counter)],...
    [Y_Min Y_Min Y_Max Y_Max Y_Min],...
    'c','FaceColor', Colors(counter,:),'FaceAlpha', 0.2,'LineStyle','none');
end
ylim([Y_Min Y_Max])
xlim([Section_Values(1) Section_Values(end)])

hold off
end