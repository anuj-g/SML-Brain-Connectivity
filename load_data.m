X_1= csvread('features_x.csv');

%X_1=featuresx();
X_2=dataset();
%X=horzcat(X_1,X_2);
X=X_2;
%X=X_1;
%{
X=data(1:70,1:19);
%}
Labels =csvread('features_y.csv');


%disp(y);
%{
num_shuffles =10;
num_folds=10;
%Fit = zeros(num_shuffles,num_folds);
%disp(size(B));
%disp(B);

[b,dev,stats] = glmfit(X,Labels,'binomial','logit');
[yfit,dylo,dyhi] = glmval(b,X,'logit',stats);
my_y=zeros(1,114);
for k =1:114
    %disp(yfit(k));
    if yfit(k)>0.5
        my_y(k)= 1;
    else
        my_y(k)=0;
    end
end
%disp(my_y);
count=0;
for k =1:114
    %disp(yfit(k));
   % if y(k)==my_y(k)
        count=count+1;
   % end
end
%disp(count/114);

sum=0;

for i = 1:1
    C = crossvalind('Kfold', Labels, 10);
    Test = (C == i); 
    Train = ~Test;                  
    SVMStruct = svmtrain ( X (Train,:), Labels (Train,:));
    [b,dev,stats] = glmfit(X(Train,:),Labels(Train,:),'binomial','logit');
    yhat=glmval(b,X(Test,:), 'logit' );
    class1=yhat<0.5;
    class2=yhat>0.5;
    disp class1
    size(yhat)
    size(yfit);
    my_y=zeros(size(yfit))';
    my_y;
    for k =1:size(yfit)
    disp(yfit(k));
        if yfit(k)>0.5
           my_y(k)= 1;
        else
            my_y(k)=0;
        end
    end

    Result = svmclassify(SVMStruct, X (Test,:));
    sum=sum+size(find( Result == Labels(Test,:)))/size(Labels(Test,:));
    
end
disp('my cross val');
disp(sum*10);
%}
disp('svm');
SVMModel = fitcsvm(X,Labels,'Standardize',true);
sum_2=0;
for i=1:10
    
    CVSVMModel = crossval(SVMModel,'KFold',10);
    FirstModel = CVSVMModel.Trained{1};
    
    sum_2= sum_2+kfoldLoss(CVSVMModel);
end
1 - (sum_2/10)

%{
disp('svm labels')
label= predict(SVMModel,X);
label
size(label)
mycount=0;
%size(label==Labels)
for i=1:44
   % disp(label{i,1});
    %disp(Labels(i,1));
    if label(i,1) == Labels(i,1);
        mycount = mycount+1;
    end
end
disp('mycount')
disp(mycount/44);
%}
sum_linear=0;


disp('linear')
for i=1:10
    LinearModel = fitclinear(X,Labels,'KFold',10);
   % y=predict(LinearModel.Trained{1},X)
    sum_linear=sum_linear+ kfoldLoss(LinearModel);
end

1- (sum_linear/10)
%{
disp('log reg');
LogRegModel = fitglm(X,Labels);
%yeval= feval(LogRegModel,X)
%CVLogRegModel = crossval(LogRegModel);
1- kfoldLoss(LogRegModel)
%}

 %{
[yfit,dylo,dyhi] = glmval(b,X,'logit',stats);
disp(yfit);
plot(x, y,'o',x,yfit,'-','LineWidth',2)
 %}
%{
for j = 1:num_shuffles
    indices = crossvalind('Kfold',Labels,num_folds);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        [b,dev,stats] = glmfit(X(train,:),Labels(train),'binomial','logit'); % Logistic regression
        Fit(j) =   glmval(b,X(test,:),'logit')';
    end
end
%}
%disp(Fit);