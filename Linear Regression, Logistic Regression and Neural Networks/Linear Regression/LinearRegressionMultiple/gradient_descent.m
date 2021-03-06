function theta = gradient_descent(X, y, theta, alpha, iterations, do_plot)
    %GRADIENT_DESCENT do Gradient Descent for a given X, y, theta, alpha
    %for a specified number of iterations

    %if less than 6 arguments was given, then set do_plot to be false
    if nargin < 6
        do_plot = false;
    end
    if(do_plot)
        plot_hypothesis(X, y, theta);
        drawnow; pause(0.1); 
    end

    m = size(X, 1); %number of training examples
    cost_vector = []; %will store the results of our cost function

    for it = 1:iterations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gradient descent
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         theta_0 = theta(1);
         theta_rest= theta(2:end); %modified

        %update theta(1) and store in temporary variable theta_0
        sigma = 0.0;
            
            for i = 1:m
             hypothesis = calculate_hypothesis(X, theta, i); 
             output = y(i);
             sigma = sigma + (hypothesis - output);
            end
        
        theta_0 = theta_0 - ((alpha * 1.0) / m) * sigma;
        
        sigma = 0.0;


         for i = 1:m
            hypothesis = calculate_hypothesis(X, theta, i); 
            output = y(i);
            sigma = sigma + (hypothesis - output) * X(i, 2:end); %modified
         end
        theta_rest = theta_rest - ((alpha * 1.0) / m) * sigma; %modified
       
         %update theta
        theta = [theta_0, theta_rest]; %modified

        %update cost_vector
        cost_vector = [cost_vector; compute_cost(X, y, theta)];

        if do_plot
            plot_hypothesis(X, y, theta);
            drawnow; pause(0.1); 
        end
    end

    disp 'Gradient descent is finished.'
        
    if do_plot
        disp 'Press enter!'
        pause;
    end

    plot_cost(cost_vector);
        
    disp 'Press enter!';
    pause;
end