

% cell1 = test_c3.Veg_type_2(3);
% str1 = string(cell1);
% result2 = strsplit(str1);
% %test_c3.Veg_type_2(3) = result2;
% 
% table2 = table(POI.Class,'VariableNames',{'Class'});
% numc = size(POI.Class,1);
% table2.Veg_type = strings(numc,1);
% for n=1:numc
% table2.Veg_type(n) = strsplit(string(test_c3.Veg_type_1(n)));
% end

%result2 = strsplit(string(test_c3.Veg_type_2(n)));

inT = readtable('VPAR_T.xlsx');
Headers= inT.Properties.VariableNames;
Veg_n = contains(Headers,'Veg_type');
Veg_n_num = sum(Veg_n);
Num_vt = size(inT,1); %Number of landcovers
Veg_cell = cell(Num_vt,Veg_n_num);
inT_Vo = inT(:,Veg_n);
for n=1:Veg_n_num
    Veg_cell(:,n) = inT_Vo{:,n};
end

%Same for crowns
Cro_n = contains(Headers,'Ccrowns');
Cro_n_num = sum(Cro_n);
IntT_cc = [inT{:,Cro_n}]; %array of values
Cro_cell = cell(Num_vt,1);
for v=1:Num_vt
    Cro_cell(v) = {IntT_cc(v,:)}; %This works
end


% Cro_cell = cell(Num_vt,Cro_n_num);
% %inT_Co = inT(:,Cro_n);
% for c=1:Cro_n_num
%      Cro_cell(:,c) = num2cell(IntT_cc(:,c));
% end

inT_2 = inT;
for v=1:Num_vt
    inT_2.Veg_type{v} = Veg_cell(v,:); %NOTE use of curly braces to take cell array into table
    inT_2.Ccrowns{v} = Cro_cell(v,:);
end
%inT_2.Veg_type = string(inT_2.Veg_type);

outT = table(inT_2.Class,inT_2.Veg_type,Cro_cell,inT_2.Cwat,inT_2.Curb,inT_2.Crock,inT_2.Cbare);
outT.Properties.VariableNames = {'Class','Veg_type','Ccrowns','Cwat','Curb','Crock','Cbare'};


string(outT.Veg_type{1}); %Works for one but not the table