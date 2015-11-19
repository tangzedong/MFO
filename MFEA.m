function data_MFEA = MFEA(Tasks,pop,gen,selection_process,rmp,p_il)
    clc    
    tic       
    data_MFEA.success = 1;
    if mod(pop,2) ~= 0
        pop = pop + 1;
    end   
    no_of_tasks=length(Tasks);
    if no_of_tasks <= 1
        error('At least 2 tasks required for MFEA');
    end
    funvalue = -1*ones(no_of_tasks,gen*pop);
    ntotalcall = ones(1,no_of_tasks);
    D=zeros(1,no_of_tasks);
    for i=1:no_of_tasks
        D(i)=Tasks(i).dims;
    end
    D_multitask=max(D);
    options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton'); 
    
    fnceval_calls = 0;
    calls_per_individual=zeros(1,pop);
    EvBestFitness = zeros(no_of_tasks,gen);
    TotalEvaluations=zeros(1,gen);
    bestobj=inf*(ones(1,no_of_tasks));
        
    for i = 1 : pop
        population(i) = Chromosome();
        population(i) = initialize(population(i),D_multitask);
        population(i).skill_factor=0;
    end
    parfor i = 1 : pop
        [population(i),calls_per_individual(i)] = evaluate(population(i),Tasks,p_il,no_of_tasks,options);
    end
    for i = 1:pop
        for t = 1:no_of_tasks
            funvalue(t,ntotalcall(t)) = population(i).factorial_costs(t);
        end
        ntotalcall = ntotalcall + 1;
    end
    fnceval_calls=fnceval_calls + sum(calls_per_individual);
    TotalEvaluations(1)=fnceval_calls; %TotalEvaluations记录目标函数的执行次数
    
    factorial_cost=zeros(1,pop);
    for i = 1:no_of_tasks
        for j = 1:pop
            factorial_cost(j)=population(j).factorial_costs(i);
        end
        [xxx,y]=sort(factorial_cost);
        population=population(y);
        for j=1:pop
            population(j).factorial_ranks(i)=j; 
        end
        bestobj(i)=population(1).factorial_costs(i);   %bestobj(i)记录目标函数(i)的最好值
        EvBestFitness(i,1)=bestobj(i);              %EvBestFitness(i)记录每一代的(i)最好值
        bestInd_data(i)=population(1);              %bestInd_data(i)记录最好目标函数值对应的个体
    end
    for i=1:pop
        [xxx,yyy]=min(population(i).factorial_ranks);
        x=find(population(i).factorial_ranks == xxx);
        equivalent_skills=length(x);
        if equivalent_skills>1
            population(i).skill_factor=x(1+round((equivalent_skills-1)*rand(1)));
            tmp=population(i).factorial_costs(population(i).skill_factor);
            population(i).factorial_costs(1:no_of_tasks)=inf;
            population(i).factorial_costs(population(i).skill_factor)=tmp;
        else
            population(i).skill_factor=yyy;
            tmp=population(i).factorial_costs(population(i).skill_factor);
            population(i).factorial_costs(1:no_of_tasks)=inf;
            population(i).factorial_costs(population(i).skill_factor)=tmp;
        end
    end
        
    mu = 10; % Index of Simulated Binary Crossover (tunable)
    sigma = 0.02; % standard deviation of Gaussian Mutation model (tunable)
    generation=0;
    eps = 0.1;
%     figure;
    while generation <= gen %bestobj(1) > eps%
        generation = generation + 1;
        indorder = randperm(pop);
        count=1;
        for i = 1:pop
            ppt(:, i) = evaluate_TEST(population(i),Tasks,p_il,no_of_tasks,options);
        end
        
%         clf;
%         plot(gca, ppt(1,:), ppt(2,:),'b*');
%         title(['Generation = ', num2str(generation)]);
%         xlabel('Task1 objective value');
%         ylabel('Task2 objective value');
        listcc = zeros(1,pop);
        nnnn = 0;
        for i = 1 : pop/2     
            p1 = indorder(i);
            p2 = indorder(i+(pop/2));
                        
            child(count)=Chromosome();
            child(count+1)=Chromosome();
            r1 = rand(1);
            if (population(p1).skill_factor == population(p2).skill_factor) || (r1<rmp)            
                u = rand(1,D_multitask);
                cf = zeros(1,D_multitask);
                cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
                cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));
                child(count) = crossover(child(count),population(p1),population(p2),cf);
                child(count+1) = crossover(child(count+1),population(p2),population(p1),cf);
% % % % % % % %                 if rand(1) < 0.1 
% % % % % % % %                     child(count)=mutate(child(count),child(count),D,sigma/2);
% % % % % % % %                     child(count+1)=mutate(child(count+1),child(count+1),D,sigma/2);
% % % % % % % %                 end        
                sf1=1+round(rand(1));
                sf2=1+round(rand(1));
                if sf1 == 1
                    child(count).skill_factor=population(p1).skill_factor;
                else
                    child(count).skill_factor=population(p2).skill_factor;
                end
                if sf2 == 1
                    child(count+1).skill_factor=population(p1).skill_factor;
                else
                    child(count+1).skill_factor=population(p2).skill_factor;
                end
%                 if (population(p1).skill_factor ~= population(p2).skill_factor) && (r1 < rmp)
%                     listcc(count) = 1;
%                     listcc(count+1) = 1;
%                     nnnn=nnnn+2;
%                     tmpy0(count, :) = evaluate_TEST(child(count), Tasks, p_il, no_of_tasks, options);
%                     tmpy0(count+1, :) = evaluate_TEST(child(count+1), Tasks, p_il, no_of_tasks, options);
%                     hold on
%                     plot(gca, tmpy0(count, 1),tmpy0(count, 2), 'ok');
%                     hold on
%                     plot(gca, tmpy0(count+1, 1),tmpy0(count+1, 2), 'ok');
%                 end
            else
                child(count)=mutate(child(count),population(p1),D_multitask,sigma);
                child(count).skill_factor=population(p1).skill_factor;
                child(count+1)=mutate(child(count+1),population(p2),D_multitask,sigma);
                child(count+1).skill_factor=population(p2).skill_factor;
            end
            count=count+2;
        end
        parfor i = 1 : pop            
            [child(i),calls_per_individual(i)] = evaluate(child(i),Tasks,p_il,no_of_tasks,options);           
        end
        for i = 1:pop
            skillgroup = child(i).skill_factor;
            funvalue(skillgroup,ntotalcall(skillgroup)) = child(i).factorial_costs(skillgroup);
            ntotalcall(skillgroup) = ntotalcall(skillgroup) + 1;
        end
%         for i = 1:pop
%             if listcc(i)
%                     tmpy = evaluate_TEST(child(i),Tasks,p_il,no_of_tasks,options);
%                     hold on
%                     plot(gca, tmpy(1), tmpy(2), 'or');
%                     hold on
%                     line([tmpy0(i,1), tmpy(1)], [tmpy0(i,2), tmpy(2)],'Color',[0.3, 0.6, 0.1], 'Parent', gca);
%             end
%         end
        fnceval_calls=fnceval_calls + sum(calls_per_individual);
        TotalEvaluations(generation)=fnceval_calls;
        
        intpopulation(1:pop)=population;
        intpopulation(pop+1:2*pop)=child;
        factorial_cost=zeros(1,2*pop);
        for i = 1:no_of_tasks
            for j = 1:2*pop
                factorial_cost(j)=intpopulation(j).factorial_costs(i);
            end
            [xxx,y]=sort(factorial_cost);
            intpopulation=intpopulation(y);
            externalpop((i-1)*10+1:(i-1)*10+10)=intpopulation(y(1:10));
            for j=1:2*pop
                intpopulation(j).factorial_ranks(i)=j;
            end
            if intpopulation(1).factorial_costs(i)<=bestobj(i)
                bestobj(i)=intpopulation(1).factorial_costs(i);
                bestInd_data(i)=intpopulation(1);
            end
            EvBestFitness(i,generation)=bestobj(i);
            data_MFEA.bestindx = y(1);
            data_MFEA.nfc = y(1) + generation * pop;
        end
        for i=1:2*pop
            [xxx,yyy]=min(intpopulation(i).factorial_ranks);
            intpopulation(i).skill_factor=yyy;
            intpopulation(i).scalar_fitness=1/xxx;
        end   
        
        
        
        if strcmp(selection_process,'elitist')
            [xxx,y]=sort(-[intpopulation.scalar_fitness]);
            intpopulation=intpopulation(y);
            population=intpopulation(1:pop);            
        elseif strcmp(selection_process,'roulette wheel') %采用赌轮盘方式选择，同时保证每个skillgroup有相同数目的个体
            for i=1:no_of_tasks
                skill_group(i).individuals=intpopulation([intpopulation.skill_factor]==i);
            end
            count=0;
            while count<pop
                count=count+1;
                skill=mod(count,no_of_tasks)+1;
                population(count)=skill_group(skill).individuals(RouletteWheelSelection([skill_group(skill).individuals.scalar_fitness]));
            end     
        end
                
%         for i = 1:length(externalpop)
%             ppt2(:, i) = evaluate_TEST(externalpop(i),Tasks,p_il,2,options);
%         end
%         ppt2n(1, :) = (ppt2(1, :) - min(ppt2(1,:)))./ (max(ppt2(1, :)) - min(ppt2(1,:)));
%         ppt2n(2, :) = (ppt2(2, :) - min(ppt2(2,:)))./ (max(ppt2(2, :)) - min(ppt2(2,:)));
%         [idxppt2, Cppt2, ~, Dppt2] = kmeans(ppt2n', 2, 'MaxIter', 10);
%         if Cppt2(1, 1) > Cppt2(1, 2) && Cppt2(2, 1) < Cppt2(2, 2)
%             tmp = idxppt2;
%             idxppt2(tmp == 1) = 2;
%             idxppt2(tmp == 2) = 1;
%         end
%         lecor = sum((Cppt2(1, :) - Cppt2(2, :)).^2) / (max(Dppt2(idxppt2==1, 1)) + max(Dppt2(idxppt2==2, 2)));
%         clf
%         plot(ppt2(1,idxppt2==1), ppt2(2,idxppt2==1),'*')
%         hold on
%         plot(ppt2(1,idxppt2==2), ppt2(2,idxppt2==2),'*')
        
        disp(['MFEA Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj)]);
    end
    if generation >= gen
        data_MFEA.success = 0;
    end
    data_MFEA.wall_clock_time=toc;
    data_MFEA.EvBestFitness=EvBestFitness;
    data_MFEA.bestInd_data=bestInd_data;
    data_MFEA.TotalEvaluations=TotalEvaluations;
    data_MFEA.funvalue = funvalue;
    save('data_MFEA','data_MFEA');
    for i=1:no_of_tasks
        figure(i)
        hold on
        plot(EvBestFitness(i,:));
        xlabel('GENERATIONS')
        ylabel(['TASK ', num2str(i), ' OBJECTIVE'])
        legend('MFEA')
    end
    tmp = clock;
    save(['mat', num2str(tmp(4)), num2str(tmp(5))], 'EvBestFitness');
end