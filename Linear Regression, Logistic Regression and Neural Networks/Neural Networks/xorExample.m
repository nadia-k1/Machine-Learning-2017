
training_set_input = [
    0,0;
    0,1;
    1,0;
    1,1
];

training_set_output = [
    0;
    1;
    1;
    0
];
%NOR logical function
%training_set_output = [
%    1;
%    0;
%    0;
%    0
%];
%AND logical function
%training_set_output = [
%    0;
%    0;
%    0;
%    1
%];

[errors,nn,training_errors,test_errors] = NeuralNetwork.train(training_set_input,training_set_output,2,10000,1);
NeuralNetwork.test_xor(training_set_input,training_set_output,nn);
figure()
plot(errors)
xlabel('iterations')
ylabel('cost')
display('press enter:')
pause