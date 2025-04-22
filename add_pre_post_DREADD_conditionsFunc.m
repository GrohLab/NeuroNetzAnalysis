function [Conditions] = add_pre_post_DREADD_conditionsFunc(Conditions)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

ITP=77000000;
Post_Injection=struct;
Pre_Injection=struct;
post=@(x)x>ITP;
% Conditions(9) = [];
for i = 3:length(Conditions)
    for d=1:length(Conditions(i).Triggers)
        if post(Conditions(i).Triggers(d,:))
            Post_Injection(i-2).name = strcat(Conditions(i).name,'J60');
            Post_Injection(i-2).Triggers(d,:) = Conditions(i).Triggers(d,:);
        else
            Pre_Injection(i-2).name = Conditions(i).name;
            Pre_Injection(i-2).Triggers(d,:) = Conditions(i).Triggers(d,:);

        end
    end
end
%%
for i = 1:length(Post_Injection)
    Post_Injection(i).Triggers = Post_Injection(i).Triggers(any(Post_Injection(i).Triggers, 2), :);
end

for i = 1:length(Pre_Injection)
    Pre_Injection(i).Triggers = Pre_Injection(i).Triggers(any(Pre_Injection(i).Triggers, 2), :);
end

Conditions=[Conditions(1:2),Pre_Injection,Post_Injection];
end