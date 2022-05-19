function hypothesis = calculate_hypothesis(X, theta, training_example)
    %CALCULATE_HYPOTHESIS This calculates the hypothesis for a given X,
    %theta and specified training example
     

 hypothesis = 0.0;
    %modified
    for t = 1:length(theta)
        hypothesis = hypothesis + theta(t) * (X(training_example,2) ^ (t - 1));
    end

end

