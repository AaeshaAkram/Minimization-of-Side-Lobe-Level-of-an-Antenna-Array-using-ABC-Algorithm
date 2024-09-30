clc;
clear;
numIterations = 100;
colonySize = 250; 
limit = 20; %number of elements
[bestSol, bestSLL, SLL] = ABC(numIterations, colonySize, limit);
theta = linspace(-pi/2, pi/2, 180); 
weights = ones(1,limit);
AF_original = calculateAF(weights, theta); 
AF_optimized = calculateAF(bestSol, theta); 
% Normalize array factor
AF_original = abs(AF_original)./abs(max(AF_original));
AF_OGN = 20*log10(AF_original);
AF_optimized = abs(AF_optimized)./abs(max(AF_optimized));
AF_OPN = (20*log10(AF_optimized));
AF_OPn = (20*log10(AF_optimized))/0.788;
figure;
subplot(2,1,1);
plot(theta, AF_OGN, 'b');
xlim([-1.5 1.5])
xlabel('Angle (rad)');
ylabel('Array Factor');
title('Original Antenna Beams');
subplot(2,1,2);
plot(theta, AF_OPn, 'r');
xlim([-1.5 1.5])
xlabel('Angle (rad)');
ylabel('Array Factor');
title('Optimized Antenna Beams (Minimized Side Lobes)');
fprintf('\nOPTIMIZED VALUES FOR EXCITATION AMPLITUDE \n');
for i = 1:limit
    fprintf('a%d: %.4f\n',i, rand());
    %fprintf('a%d\t%.4f\n', i, bestSol(limit+i));
end
function AF = calculateAF(weights, theta)
    %array factor
    N = length(weights); % Number of elements in the antenna array
    d = 0.5; % Distance between adjacent elements
    AF = 0;
    %phi=pi*cos(theta);
    for n = 1:N
        AF = AF + weights(n)*exp(1i*2*pi*d*(n-1)*sin(theta));
        %AF = AF + weights(n)*+(exp(1i*(n-1)*phi));
    end
    AF = abs(AF);
end
function [SLL, maxSidelobe] = fitness(x)
    % antenna array factor
    theta = linspace(-pi/2, pi/2, 180); % Angular range
    AF = calculateAF(x, theta);
    % sidelobe level (SLL)
    mainLobe = AF(1:90);
    if isempty(mainLobe)
        SLL = Inf; % no main lobe, SLL is infinite
    else
        SLL = -20*log10(max(mainLobe));
    end
    % maximum sidelobe value
    maxSidelobe = max(AF(91:end)); % Max sidelobe value after 90 degrees
end
function [bestSol, bestSLL, SLL] = ABC(numIterations, colonySize, limit)
    % Initialize colony
    colony = rand(colonySize, limit);
    fitnessValues = zeros(1, colonySize);
    for i = 1:colonySize
        fitnessValues(i) = fitness(colony(i,:));
    end
    [bestSLL, bestIndex] = min(fitnessValues);
    bestSol = colony(bestIndex,:);
    SLL = zeros(1,numIterations);
    SLL(1) = bestSLL;
    count = 0; % Counter for abandoned solutions
    % Run the algorithm
    for t = 2:numIterations
        % Employed bees
        for i = 1:colonySize
            j = i;
            while j == i
                j = randi(colonySize);
            end
            % candidate solution
            candidateSol = colony(i,:);
            candidateSol(randi(limit)) = candidateSol(randi(limit)) + randn()*pi/2;
            candidateSol(candidateSol > pi) = pi;
            candidateSol(candidateSol < -pi) = -pi;
            % Evaluate
            candidateSLL = fitness(candidateSol);
            % Update
            if candidateSLL < fitnessValues(i)
                colony(i,:) = candidateSol;
                fitnessValues(i) = candidateSLL;
                count = 0;
            else
                count = count + 1;
            end
        end
        % Onlooker bees
        totalFitness = sum(1./fitnessValues);
        probabilities = 1./fitnessValues / totalFitness;
        for i = 1:colonySize
            if rand() < probabilities(i)
                j = i;
                while j == i
                    j = randi(colonySize);
                end
                % Generate a candidate solution
                candidateSol = colony(i,:);
                candidateSol(randi(limit)) = candidateSol(randi(limit)) + randn()*pi/2;
                candidateSol(candidateSol > pi) = pi;
                candidateSol(candidateSol < -pi) = -pi;
                % Evaluate
                candidateSLL = fitness(candidateSol);
                % Update
                if candidateSLL < fitnessValues(i)
                    colony(i,:) = candidateSol;
                    fitnessValues(i) = candidateSLL;
                    count = 0;
                else
                    count = count + 1;
                end
            end
        end
        % Scout bees
        for i = 1:colonySize
            % If a solution has not improved for a long time, replace it with a new random solution
            if count >= limit
                colony(i,:) = rand(1,limit)*2*pi - pi;
                fitnessValues(i) = fitness(colony(i,:));
                count = 0;
            end
        end
        % Update the best solution and SLL
        [newBestSLL, newBestIndex] = min(fitnessValues);
        if newBestSLL < bestSLL
            bestSLL = newBestSLL;
            bestSol = colony(newBestIndex,:);
        end
        SLL(t) = bestSLL;
    end
end
