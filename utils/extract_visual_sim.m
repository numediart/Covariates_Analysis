function [visualSimilarity] = extract_visual_sim(PATH_TO_SIMILARIITY, PATH_TO_ITEMS)

%% Add visual similarity values to trialinfo
% similarity_table = readtable('C:\Users\luca-\OneDrive - UMONS\_PhD\ARC\visual_similarity.csv');
% items = readtable('C:\Users\luca-\OneDrive - UMONS\_PhD\ARC\cible_amorce.csv');
similarity_table = readtable(PATH_TO_SIMILARIITY);
items = readtable(PATH_TO_ITEMS);
similarity_table = similarity_table(1:19,:);
items = items(1:19,:);
visualSimilarity = items(:,1:3);
visualSimilarity.Properties.VariableNames = {'target','primer','similarity'};

itemMapping = [7 8 5 6 3 4 1 2 11 12 9 10];

tmpTable = similarity_table;
for i=2:2:12
    for j = 1:19
        tmp = strrep(similarity_table.(['Var' num2str(i)]){j},',','.');
%         similarity_table.(['Var' num2str(i)]){j} = str2double(tmp);
        tmpTable.(['Var' num2str(i)]){j} = str2double(tmp);
        similarity_table.(['Var' num2str(i-1)]){j} = str2double(regexp(tmpTable.(['Var' num2str(i-1)]){j},'\d*','Match'));
    end
    [~,idx] = sort(cell2mat(similarity_table.(['Var' num2str(i-1)])));
    similarity_table.(['Var' num2str(i-1)]) = similarity_table.(['Var' num2str(i-1)])(idx);
    similarity_table.(['Var' num2str(i)]) = similarity_table.(['Var' num2str(i)])(idx);
end
clear tmp

for i=1:2:11
    for j = 1:19
%         tmp = regexp(similarity_table.(['Var' num2str(i)]){j},'\d*','Match');
%         tmp = similarity_table.(['Var' num2str(i)]){j};
%         idx = str2double(tmp{1});
        idx=j;
        visualSimilarity{(ceil(i/2)-1)*19+idx,1:2} = items{idx,itemMapping(i):itemMapping(i+1)};
        visualSimilarity{(ceil(i/2)-1)*19+idx,3} = {str2double(strrep(similarity_table{idx,i+1},',','.'))};
    end
end