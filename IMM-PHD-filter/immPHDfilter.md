# IMM-PHD Filter Implementation Ideas

Due to my work commitments, I haven't had time to work on this project as much as I'd like. Below are some ideas for implementing the IMM-PHD filter. If there are better solutions, please let me know.

## Motion Model Selection

1. **Constant Velocity (CV) Model**: Look up to V.O's codes in gen_model.
2. **Constant Acceleration (CA) Model**: This model accounts for acceleration, which can be significant in highly maneuverable targets.
3. **Coordinated Turn (CT) Model**: Useful for targets executing turns.

For each model, we need to define the state transition matrix (F), process noise (Q), and potentially the control input matrix (B) if external inputs are considered.

### Constant Acceleration Model

```matlab
model.F_CA = [1 model.T (model.T^2)/2; 0 1 model.T; 0 0 1];
model.F_CA = blkdiag(model.F_CA, model.F_CA);  % Assuming 2D motion
model.sigma_a = 2;  % Example value, adjust based on scenario
model.Q_CA = model.sigma_a^2 * [model.T^4/4 model.T^3/2 model.T^2/2; model.T^3/2 model.T^2 model.T; model.T^2/2 model.T 1];
model.Q_CA = blkdiag(model.Q_CA, model.Q_CA);  % Assuming 2D motion

### Coordinated Turn Model

```matlab
model.F_CT = [1 sin(model.T*omega)/omega 0 -(1-cos(model.T*omega))/omega; 0 cos(model.T*omega) 0 -sin(model.T*omega); 0 (1-cos(model.T*omega))/omega 1 sin(model.T*omega)/omega; 0 sin(model.T*omega) 0 cos(model.T*omega)];
model.omega = 0.1;  % Turn rate, rad/s
model.sigma_turn = 0.1;
model.Q_CT = model.sigma_turn^2 * eye(4);  % Simplified Q
```

# Initialize the Filter Parameters

```matlab
% Initialize IMM probabilities
model.mu = [1/3, 1/3, 1/3];  % Uniform prior

% Transition matrix
model.pi = [0.7, 0.15, 0.15;
            0.15, 0.7, 0.15;
            0.15, 0.15, 0.7];
```

# Subfunctions

## Ground Truth Generation Module

```matlab
function truth = gen_truth_with_maneuvers(model)
    % Define maneuver intervals and models
    maneuvers = [1, 20, 40; 1, 2, 3];  % Start time and model index (1-CV, 2-CA, 3-CT)
    current_model = model.F;  % Start with CV
    current_Q = model.Q;

    for k = 1:truth.K
        if any(maneuvers(1,:) == k)
            model_idx = maneuvers(2, maneuvers(1,:) == k);
            switch model_idx
                case 1
                    current_model = model.F;
                    current_Q = model.Q;
                case 2
                    current_model = model.F_CA;
                    current_Q = model.Q_CA;
                case 3
                    current_model = model.F_CT;
                    current_Q = model.Q_CT;
            end
        end

        % Generate state based on current model
        % Existing code to propagate state... refer to V.O's
    end
end
```

## Mixing Components

Before the prediction step, mix the components based on the transition probabilities. This mixes the Gaussian components across different motion models to account for the uncertainty in motion model states.

```matlab
function [mixed_w, mixed_m, mixed_P] = mix_components(w, m, P, model)
    num_models = length(model.mu);
    num_components = length(w);
    mixed_w = zeros(num_components, num_models);
    mixed_m = zeros(size(m, 1), num_components, num_models);
    mixed_P = zeros(size(P, 1), size(P, 2), num_components, num_models);
    
    for j = 1:num_models
        for i = 1:num_models
            mixing_prob = model.pi(i, j) * model.mu(i);
            idx = (j - 1) * num_components + 1 : j * num_components;
            mixed_w(idx, j) = mixing_prob * w;
            mixed_m(:, idx, j) = m;
            mixed_P(:, :, idx, j) = P;
        end
    end
    mixed_w = sum(mixed_w, 2);
    mixed_m = reshape(mixed_m, size(m, 1), []);
    mixed_P = reshape(mixed_P, size(P, 1), size(P, 2), []);
end
```

## Prediction and Update

Modify the prediction and update functions to handle multiple models. During prediction, use the appropriate model matrices (F and Q) for each model.

```matlab
function [w_pred, m_pred, P_pred] = predict_phd(w, m, P, model)
    [w_mix, m_mix, P_mix] = mix_components(w, m, P, model);
    for j = 1:size(m_mix, 2)
        [m_pred(:, j), P_pred(:, :, j)] = kalman_predict_single(model.F{j}, model.Q{j}, m_mix(:, j), P_mix(:, :, j));
        w_pred(j) = model.P_S * w_mix(j);  % Survival probability
    end
end

function [w_update, m_update, P_update] = update_phd(z, model, w_pred, m_pred, P_pred)
    if isempty(z)
        % Handle misdetection
        w_update = model.Q_D * w_pred;
        m_update = m_pred;
        P_update = P_pred;
    else
        % Normal update with measurements
        % Insert normal PHD update code here
    end
end
```

## IMM Interaction

After updating, calculate the model probabilities based on the likelihood of each model given the current measurements. Then, re-normalize the model probabilities.

```matlab
function mu_new = update_model_probabilities(z, model, m_pred, P_pred)
    mu_new = zeros(size(model.mu));
    for j = 1:length(model.mu)
        likelihood = calculate_likelihood(z, m_pred(:, j), P_pred(:, :, j), model.H, model.R);
        mu_new(j) = likelihood * model.mu(j);
    end
    mu_new = mu_new / sum(mu_new);  % Normalize
end
```

## Backbone

The backbone includes the core steps of mixing, predicting, updating, and managing the Gaussian mixtures. Note that the specifics of the subfunctions like mixing, predicting, and updating will need to fit the IMM framework, which means handling multiple models.

```matlab
function est = run_imm_phd_filter(model, meas)
    % Output variables
    est.X = cell(meas.K, 1);
    est.N = zeros(meas.K, 1);
    est.L = cell(meas.K, 1);

    % Filter parameters
    filter.L_max = 100;                  % Limit on number of Gaussians
    filter.elim_threshold = 1e-5;        % Pruning threshold
    filter.merge_threshold = 4;          % Merging threshold
    filter.P_G = 0.999;                  % Gate size in percentage
    filter.gamma = chi2inv(filter.P_G, model.z_dim); % inv chi^2 dn gamma value
    filter.gate_flag = 1;                % Gating on or off 1/0
    filter.run_flag = 'disp';            % 'disp' or 'silence' for on the fly output

    est.filter = filter;

    % Initialize prior probabilities for each model
    w_update = eps * ones(model.L_birth, 1);
    m_update = repmat([0.1; 0; 0.1; 0], 1, model.L_birth);
    P_update = repmat(diag([1 1 1 1]).^2, [1, 1, model.L_birth]);
    L_update = model.L_birth;

    % Recursive filtering
    for k = 1:meas.K

        % Prediction
        [w_predict, m_predict, P_predict] = predict_imm_phd(w_update, m_update, P_update, model);

        % Append birth components
        w_predict = [model.w_birth; w_predict];
        m_predict = [model.m_birth, m_predict];
        P_predict = cat(3, model.P_birth, P_predict);

        % Gating (if enabled)
        if filter.gate_flag
            meas.Z{k} = gate_meas_gms(meas.Z{k}, filter.gamma, model, m_predict, P_predict);
        end

        % Update
        if isempty(meas.Z{k})
            w_update = model.Q_D * w_predict;
            m_update = m_predict;
            P_update = P_predict;
        else
            [w_update, m_update, P_update] = update_imm_phd(meas.Z{k}, model, w_predict, m_predict, P_predict);
        end

        % Mixture management: pruning, merging, capping
        [w_update, m_update, P_update] = gaus_prune(w_update, m_update, P_update, filter.elim_threshold);
        [w_update, m_update, P_update] = gaus_merge(w_update, m_update, P_update, filter.merge_threshold);
        [w_update, m_update, P_update] = gaus_cap(w_update, m_update, P_update, filter.L_max);

        % State extraction
        idx = find(w_update > 0.5);
        for j = 1:length(idx)
            repeat_num_targets = round(w_update(idx(j)));
            est.X{k} = [est.X{k}, repmat(m_update(:, idx(j)), 1, repeat_num_targets)];
            est.N(k) = est.N(k) + repeat_num_targets;
        end

        % Display diagnostics
        if ~strcmp(filter.run_flag, 'silence')
            disp([' time= ', num2str(k), ...
                ' #est mean=', num2str(sum(w_update), 4), ...
                ' #est card=', num2str(est.N(k), 4), ...
                ' #gaus orig=', num2str(length(w_predict)), ...
                ' #gaus elim=', num2str(length(w_update))]);
        end
    end
end
```






