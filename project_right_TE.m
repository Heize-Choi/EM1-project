%1. Impoting Data(unit=minute)
sleepTbl = readtable('Example_Sleep_Wake_Schedule_3.csv');
lightTbl = readtable('illuminance_2.csv');
sleepTbl.sleep_start = datetime(sleepTbl.sleep_start, 'InputFormat','yyyy-MM-dd HH:mm:ss');
sleepTbl.sleep_end   = datetime(sleepTbl.sleep_end,   'InputFormat','yyyy-MM-dd HH:mm:ss');
lightTbl.timestamp   = datetime(lightTbl.timestamp,   'InputFormat','yyyy-MM-dd HH:mm:ss');

startTime = dateshift(min(lightTbl.timestamp),'start','day');
minuteVec = (startTime:minutes(1):startTime+days(6)-minutes(1))'; %% 6 Days Data(Cumulative cycle 5 days + TE measurement 1 day)

%%Illumination(lux)
illumVec = interp1( ...
  datetime(lightTbl.timestamp), lightTbl.illuminance, ...
  datetime(minuteVec), 'linear','extrap' );

%%where sleep or wake
isAwake = true(size(minuteVec));
for i=1:height(sleepTbl)
  idx = minuteVec>=sleepTbl.sleep_start(i) & minuteVec<sleepTbl.sleep_end(i);
  isAwake(idx) = false;
end


%% 2. Initial valuse & Parameters
%24 hours per minute 
n_day = 1440;
t_last  = (0:n_day-1).';
tt_hours = t_last / 60;

t0 =   0;    
T  = 1440*6;   %total dat length
h  = 1;   
t = t0:h:T; 

%for H1 rythm
n0   = 0;    
xc0  = 1;   
x0   = 0;    
H10  = 0;    

beta     = 0.0075;
kappa    = 12/pi;
f_val    = 0.99729;   
tau_c    = 24.1;
K        = 0.5; 
zeta     = 0.54;
M        = 0.019513;
eta      = 0.04;
beta_Ip  = 7.83e-4;
m        = 7;
phi_on = 6.113;
phi_off = 4.352;
a = 1.0442e-3;
delta = 600;
r = 15.36;
G = 37;
b = 0.4;
gamma = 0.13;

%for TE
alpha   = 0.1;   
I0       = 9500;   
p        = 0.5;   
I1       = 100;
A = 100;
B = 10;
f_val2 = 0.0026564;
K_TE = 0.5;
Rc = 2880;
Rt(1) = 1440; 
I_max = 5;
i_val = 0.04;
as = 0.55;    
Rt = zeros(size(t));
t_a = zeros(size(t));

%% 3. functions
% for t_a
wakeIdx = find(diff(isAwake)==1) + 1;  
for k = 1:numel(wakeIdx)
    i0 = wakeIdx(k);
    iEnd = find(t >= t(i0) + 120, 1) - 1;
    if isempty(iEnd)
        iEnd = numel(t);
    end
    idx = i0:iEnd;
    t_a(idx) = t(idx) - t(i0);
end

t_a_last = t_a(end-n_day+1:end);

%%Rt
for i = 1:length(t)-1
    if isAwake(i)
        Rt(i+1) = Rt(i) - K_TE;
    else
        Rt(i+1) = Rt(i) + f_val2*(Rc - Rt(i));
    end
end

Rt_last = Rt(end-n_day+1:end); % t_last와 같은 길이
I_fun = @(tt) interp1( ...
           minuteVec, illumVec, ...
           startTime + minutes(tt), ...
           'linear', 'extrap' );

awake_fun = @(tt) logical(interp1( ...
                 minuteVec, double(isAwake), ...
                 startTime + minutes(tt), ...
                 'nearest','extrap'));

alpha_fun = @(n, I) alpha*(I./I0).^p .* (I./(I + I1));

B_fun     = @(n,x,xc,I) G.*(1 - b.*x).*(1 - b.*xc).*alpha_fun(1 - n, I);  

S_fun     = @(H1,B) double( (1  - m.*B)>=0 | (H1>=0.001) );

A_phi = @(x,xc) ...
  (a*(1 - exp(-delta * mod(phi_on  - mod(atan(x./xc),2*pi),2*pi))) ...
     ./ (1 - exp(-delta * mod(phi_on - phi_off,2*pi)))) ...
    .* ((mod(atan(x./xc),2*pi)>phi_off)|(mod(atan(x./xc),2*pi)<phi_on)) ...
  + (a*exp(-r * mod(phi_on - phi_off,2*pi))) ...
    .* (~((mod(atan(x./xc),2*pi)>phi_off)|(mod(atan(x./xc),2*pi)<phi_on)));





%% 4. Runge-Kutta            
N = numel(t);
y = zeros(4, N);
y(:,1) = [ n0; x0; xc0; H10 ];

f = @(tt, yy)[     
  % 1)n
  alpha_fun(yy(1), I_fun(tt))*(1-yy(1)) - beta*yy(1);

  % 2) ẋ
  1/kappa*( ...
     yy(3) ...
   + gamma*( yy(2)/3 + 4*yy(2).^3/3 - 256*yy(2).^7/105 ) ...
   + B_fun(yy(1),yy(2),yy(3),I_fun(tt)) ...
   - eta*M ...
  );
    % 3) ẋc
  1/kappa*( ...
     B_fun(yy(1),yy(2),yy(3),I_fun(tt))*yy(3)/3 ...
   - yy(2)*(24/(f_val*tau_c))^2 ...
   + K*B_fun(yy(1),yy(2),yy(3),I_fun(tt)) ...
   - zeta*M ...
  );
  % 4) H1̇
  -beta_Ip*yy(4) ...
  + A_phi(yy(2),yy(3)) ...
    .*(1 - m*B_fun(yy(1),yy(2),yy(3),I_fun(tt))) ...
    .* S_fun(yy(4), B_fun(yy(1),yy(2),yy(3),I_fun(tt)))
];

% RK4 loop
for i = 1:N-1
  k1   =  f( t(i),        y(:,i) );
  k2   =  f( t(i)+h/2,    y(:,i)+k1*h/2 );
  k3   = f( t(i)+h/2,    y(:,i)+k2*h/2 );
  k4   =  f( t(i)+h,      y(:,i)+k3*h   );
  y(:,i+1) = y(:,i) + (k1 + 2*k2 + 2*k3 + k4)*h/6;
end

H1 = y(4, :);
H1_last = H1(end-n_day+1 :end); %H1 rythm for day

%% 5. H1 transform for TE
w1 = 2*pi./24;
w2 = 4*pi./24;
X = [ ones(n_day,1), ...
      cos(w1*tt_hours), ...
      sin(w1*tt_hours), ...
      cos(w2*tt_hours), ...
      sin(w2*tt_hours) ];

b_coef = X \ H1_last.'; 

a0   = b_coef(1);
A1   = b_coef(2);  B1 = b_coef(3);
a1   = sqrt(A1^2 + B1^2);
phi1 = atan2(B1, A1) * 24/(2*pi);

A2   = b_coef(4);  B2 = b_coef(5);
a2   = sqrt(A2^2 + B2^2);
phi2 = atan2(B2, A2) * 24/(4*pi);

C1 = a1;
C2 = a2/a1;
P = phi1;
pp = phi2-phi1;

%% 6. TE funtion

H_fun = @(t) C1.*(cos(2*pi*((t/60)-P)/24)+(C2*cos(4*pi*((t/60)-P+pp)/24)));

C_fun = @(t) H_fun(t)./ (C1*(1+C2));


I_func = @(t) (-I_max .* ...
              exp(- (i_val .* t_a_last(t)) ./ ...
                      (-as .* C_fun(t) + f_val2 .* (Rc - Rt_last(t))))) ...
              .* (t_a_last(t) > 0);

TE = @(t) 80.* (Rt_last(t)./ Rc) + B + H_fun(t) + I_func(t);

ts = 1:1440;
TE_values = TE(ts);
TE_values(TE_values < 0) = 0;

for k = 1:length(ts)
    fprintf('TE(%4d) = %g\n', k, TE_values(k));
end

%% (add) Plots

figure;
subplot(2,1,1);
plot(minuteVec, illumVec); ylabel('illuminance_2.csv'); title('Illuminance Pattern(6 days)');
subplot(2,1,2);
plot(t, H1); xlabel('Time (min)'); ylabel('H_1'); title('H1 (6 days)');

%근사 확인
figure;
plot(tt_hours, H1_last, 'b-'); hold on;
plot(tt_hours, a0 + a1*cos(2*pi*(tt_hours-phi1)./24) + a2*cos(4*pi*(tt_hours-phi2)./24), 'r--');
xlabel('Time (h)'); ylabel('H_1'); legend('Simulated','Cosine Fit');
title('H1 model');
figure;
plot(ts, TE_values, 'LineWidth', 1.5);
xlabel('t (minutes)');
ylabel('TE(t)');
title('Task Effectiveness');
grid on;


% 전체 구간 시간 벡터 (시간 단위)
tt_hours_all = t / 60; % t는 전체 분 단위 시간 벡터

% 전체 구간 H1
H1_all = H1; % 또는 y(4,:)

% 전체 구간 코사인 근사 (전체 구간에 대해 피팅한 a0, a1, phi1, a2, phi2 사용)
cosine_fit_all = a0 + a1*cos(2*pi*(tt_hours_all-phi1)./24) + a2*cos(4*pi*(tt_hours_all-phi2)./24);

% R² 계산
ss_res_all = sum((H1_all - cosine_fit_all).^2);
ss_tot_all = sum((H1_all - mean(H1_all)).^2);
R2_all = 1 - (ss_res_all / ss_tot_all);

% 출력
fprintf('전체 구간 대표 R²: %.4f\n', R2_all);

%%
T = readtable('timelogs.csv', 'ReadVariableNames', false); % 헤더 없는 파일
test_time = T{:,1};
timestamp_raw = string(T{:,4});
error_times=T{:,2};

timestamp = datetime(timestamp_raw, 'InputFormat', 'yyyyMMdd_HHmmss');

test_edd = 300 ./ (test_time+error_times*(0.1));
te_time = startTime + minutes(ts);

figure;
yyaxis left;
plot(te_time, TE_values, 'b-', 'LineWidth', 1.5);
ylabel('Task Effectiveness (TE(t))');
ylim([0, max(TE_values)+10]);

yyaxis right;
plot(timestamp, test_edd, 'r--o', 'LineWidth', 1.2);
ylabel('300 / test time');
ylim([0, max(test_edd)+10]);

xlabel('Time');
title('TE(t)와 1 / 조도 비교');
legend('TE(t)', '300 / time', 'Location', 'best');
grid on;
%% 
t_min = minutes(timestamp - startTime);  % test_edd와 같은 길이
% 초기 추정값 [P, C2, pp, A, B, C1]

fprintf('P   = %.4f\n', P);
fprintf('C2  = %.4f\n', C2);
fprintf('pp  = %.4f\n', pp);
fprintf('A   = %.4f\n', A);
fprintf('B   = %.4f\n', B);
fprintf('C1  = %.4f\n', C1);
fprintf('I_max  = %.4f\n', I_max);
params0 = [P, C2, pp, A, B, C1, I_max];
lb = [-Inf, -Inf, -Inf, 0, -Inf, -Inf, 0];
% 모델 함수 정의 (t 대신 t_last 사용)
TE_model = @(params, t_min) ...
    params(4) * (interp1(t_last, Rt_last, t_min, 'nearest', 'extrap') / Rc) ... % A 계수
    + params(5) ... % B 상수
    + params(6) * (cos(2*pi*(t_min - params(1))/24)) ...     % C1 * cos(2pi(t - P))
    + params(2) * cos(4*pi*(t_min - params(1) + params(3))/24) ... % + C1 * C2 * cos(4pi(t - (P + pp)))
    - params(7)* exp(-((i_val * interp1(t_last, t_a_last, t_min, 'nearest', 'extrap')) ./ ...
       (-as * params(6) * (cos(2*pi*(t_min - params(1))/24) + ...
       params(2) * cos(4*pi*(t_min - params(1) + params(3))/24)) ...
       + f_val2 * (Rc - interp1(t_last, Rt_last, t_min, 'nearest', 'extrap')))));
% 초기 추정값 [P, C2, pp]
      % 혹은 [17.2, 0.5, 1.8] 등
options = optimoptions('lsqcurvefit', 'MaxFunctionEvaluations', 1000, 'MaxIterations', 500, 'Display', 'iter');  % 예시로 최대 함수 평가 횟수를 1000으로 설정

% 피팅 수행
params_fit = lsqcurvefit(TE_model, params0, t_min, test_edd, lb, [], options);
disp(params_fit); %params_fit=[a1, a2, phi1, phi2, A, B]
% 모델 근사 결과
fitted_TE = TE_model(params_fit, t_min);

% 시각화
figure;
plot(t_min, test_edd, 'ko'); hold on;
plot(t_min, fitted_TE, 'r-', 'LineWidth', 1.5);
xlabel('Time (minutes)'); ylabel('TE');
legend('실측 TE', '모델 근사');
title('모델 피팅 결과');
grid on;

%%
P  = params_fit(1);
C2 = params_fit(2);
pp = params_fit(3);
A  = params_fit(4);
B  = params_fit(5);
C1 = params_fit(6);
I_max = params_fit(7);

t_min = 1:1440;  % 7일 동안 (10080분)
TE = @(t) A.* (Rt_last(t)./ Rc) + B + H_fun(t) + I_func(t);
% TE 값을 계산
fixed_TE_values = TE(t_min);
fixed_TE_values(fixed_TE_values < 0) = 0;

% 결과 시각화
figure;
yyaxis left;
plot(te_time, fixed_TE_values, 'b-', 'LineWidth', 1.5);
ylabel('Task Effectiveness (TE(t))');
ylim([0, max(TE_values)+10]);

yyaxis right;
plot(timestamp, test_edd, 'r--o', 'LineWidth', 1.2);
ylabel('300 / test time');
ylim([0, max(TE_values)+10]);

xlabel('t (minutes)');
title('Task Effectiveness');
grid on;
%%
mse = mean((test_edd - fitted_TE).^2);
fprintf('MSE (Mean Squared Error) = %.4f\n', mse);
