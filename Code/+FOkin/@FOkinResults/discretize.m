function disc = discretize(self)

assert((self.level == 2), ['Discretization is possible only on the ',...
    'results of ''discretize'' operation.'])

    function setval
        % nested function for calculation of centres
        values = xii(feature);
        abs_values = abs(values);
        value = sum(values);
        abs_value = sum(abs_values);
        col = sum(feature.*abs_values)/abs_value;
        col = round(col);
        centers(col, ii) = value;
    end

x = self.x;
opt = self.options;
mask = (abs(x) > opt.nonzero_limit);
x = x.*mask;
param = self.group_param;

if all(all(x == 0))
    threads = [];
    average_tau = [];
else
    tau = self.tau;
    if ~isempty(self.group_objval)
        % The result is based on a group-LASSO optimization.
        % If an element of the group has nonzero value of x at a
        % point then all other elements have also nonzero value of
        % x at the same point. Hence the discretization can be
        % based on the nenzero values of abs(x) summed over the
        % groups.
        glob_x = sum(abs(x), 2);
        pos = find(glob_x > 0); % positions of the nonzero elements
        feature_lims = [];
        idx = 1;
        stop = false;
        k = length(pos);
        while ~stop
            lims.first = pos(idx);
            old = lims.first;
            new = old;
            while ((new - old) <= opt.min_feature_gap)
                old = new;
                idx = idx + 1;
                if idx > k
                    stop = true;
                    break
                end
                new = pos(idx);
            end
            lims.last = old;
            feature_lims = [feature_lims, lims]; %#ok<AGROW> 
        end
        
        k = length(feature_lims);
        for ii = k:-1:1 % backwards for preallocation
            first = feature_lims(ii).first;
            last = feature_lims(ii).last;
            feature_x = (x(first:last,:))';
            abs_feature_x = abs(feature_x);
            feature_tau = (tau(first:last))';
            val = sum(feature_x, 2);
            abs_val = sum(abs_feature_x, 2);
            thread = FOkin.Thread;
            thread.param = param;
            thread.tau = ((abs_feature_x*feature_tau)./abs_val)';
            thread.val = (val)';
            threads(ii) = thread;
        end
    else
        % Calculate first the center and the summed value of the features.
        [n, p] = size(x);
        centers = zeros(n, p);
        for ii=1:p
            xii = x(:,ii);
            pos = find(abs(xii) > 0);
            k = length(pos);
            if (k > 0)
                feature = pos(1);
                if (k == 1)
                    centers(feature, ii) = xii(feature);
                end
                for jj=2:k
                    step = pos(jj) - feature(end);
                    test = (step > opt.min_feature_gap + 1);
                    if test
                        setval
                        feature = pos(jj);
                    else
                        feature = [feature; pos(jj)]; %#ok<AGROW>
                    end
                    if (jj == k )
                        setval
                    end
                end
            end
        end
        
        % Sort the centers into threads along the values of group_param.        
        [nrows, ncols] = size(centers);
        centers_orig = centers;
        
        threads = [];
        fragments = [];
        while 1
            test = sum(sum(abs(centers)));
            if test == 0
                break
            end
            maxval = max(max(abs(centers)));
            [row, col] = find((abs(centers) == maxval), 1);
            obj.r = row;
            obj.c = col;
            thread = obj;
            centers(row,col) = 0;            
            row_start = row;
            col_start = col;
            for step = [1, -1] % sep forward and backward from the maximum
                while 1
                    col = col + step;
                    if ((col>0) && (col <= ncols))
                        found = 0;
                        for ii = 0:opt.thread_row_region
                            if ((row+ii) <= nrows)
                                if (abs(centers(row+ii, col)) > 0)
                                    newrow = row + ii;
                                    found = 1;
                                end
                            end
                            if ((row-ii) > 0)
                                if (abs(centers(row-ii, col)) > 0)
                                    if ~found
                                        newrow = row - ii;
                                        found = 1;
                                    else
                                        % Check if this is a better point 
                                        % than that found at row + ii.
                                        newrow2 = row - ii;
                                        if (abs(centers(newrow2, col)) >...
                                                abs(centers(newrow, col)))
                                            newrow = newrow2;
                                        end
                                    end
                                end
                            end
                            if found
                                obj.r = newrow;
                                obj.c = col;
                                thread = [thread, obj]; %#ok<AGROW>
                                centers(newrow,col) = 0;
                                row = newrow;
                                break
                            end
                        end
                    else
                        break
                    end
                end
                row = row_start;
                col = col_start;
            end
            
            k = length(thread);
            rows = zeros(k,1);
            cols = zeros(k,1);
            val = zeros(1,k);
            for ii= 1:k
                obj = thread(ii);
                rows(ii) = obj.r;
                cols(ii) = obj.c;
                val(ii) = centers_orig(obj.r, obj.c);
            end
            component = FOkin.Thread;
            component.param = param(cols);
            component.tau = tau(rows);
            component.val = val;
            [component.param, idx] = sort(component.param);
            component.tau = component.tau(idx);
            component.val = component.val(idx);

            if (length(thread) >= opt.min_thread_length)
                if isempty(threads)
                    threads = component;
                else
                    threads = [threads, component]; %#ok<AGROW>
                end
            else
                if isempty(fragments)
                    fragments = component;
                else
                    fragments = [fragments, component]; %#ok<AGROW>
                end
            end
        end
        
        nt = length(threads);
        nf = length(fragments);
        if ((nt == 0) && (nf ~= 0))
            error(['Only fragment threads were found with finite tau ',...
                'value. Decrease the ''min_thread_length'' property ',...
                'in the ''options'' property of the ''FOkin'' object.'])
        end
        
        for ii=1:nf
            fragment = fragments(ii);
            % Handle the elements of the fragments individually.
            for jj = 1:length(fragment.param)
                obj.param = fragment.param(jj);
                obj.tau = fragment.tau(jj);
                obj.val = fragment.val(jj);
                for k = nt:-1:1 % bacwards for preallocation
                    thread = threads(k);
                    logtau = log10(thread.tau);
                    meantau = mean(logtau);
                    maxdev = max(abs(logtau - meantau));
                    objlogtau = log10(obj.tau);
                    if (abs(meantau - objlogtau) <= maxdev)
                        % Incorpotate the object to the thread.
                        [lo, com, hi] =...
                            overlap(thread.param, obj.param);
                        idx1 = find(com);
                        if ~isempty(idx1)
                            if (obj.param == thread.param(idx1))
                                channels(1) = find(tau ==...
                                    thread.tau(idx1));
                                channels(2) = find(tau == obj.tau);
                                vals(1) = thread.val(idx1);
                                vals(2) = obj.val;
                                val = sum(vals);
                                pos = sum(channels.*abs(vals))...
                                    /sum(abs(vals));
                                pos = round(pos);
                                obj.tau = tau(pos);
                                obj.val = val;
                            end
                        end
                        thread.param = [thread.param(lo), obj.param,...
                            thread.param(hi)];
                        thread.tau = [thread.tau(lo), obj.tau,...
                            thread.tau(hi)];
                        thread.val = [thread.val(lo), obj.val,...
                            thread.val(hi)];
                    end
                    threads(k) = thread;
                end
            end
        end
    end
end

nt = length(threads);

% Represent x0 also as a thread if any of its components is above the
% nonzero limit.
x0 = self.x0;
mask = (abs(x0) > opt.nonzero_limit);
x0= x0.*mask;
if any(x0)
    x0_thread = FOkin.Thread;
    x0_thread.param = param;
    x0_thread.tau = Inf;
    x0_thread.val = x0;
    nt = nt + 1;
    threads = [threads, x0_thread];
end

if nt > 0
    % Calculate the maximum absolute value of each thread.
    maxval = zeros(nt,1);
    for ii=1:nt
        maxval(ii) = max(abs(threads(ii).val));
    end
    
    % Sort in descending order and normalize.
    [maxval, idx] = sort(maxval,'descend' );
    abs_max_amplitude = maxval(1);
    rel_amplitude = maxval/abs_max_amplitude;
    
    % Sort the threads.
    threads = threads(idx);
    
    % Calculate the average of taus in log scale.
    average_tau = zeros(nt, 1);
    for ii = 1:nt
        logtau = log10(threads(ii).tau);
        average_tau(ii) = 10^mean(logtau);
    end
else
    abs_max_amplitude = [];
    rel_amplitude = [];
    average_tau = [];
end

disc = FOkin.FOkinDiscretized(self);
disc.discrete = FOkin.Discrete;
disc.discrete.threads = threads;
disc.discrete.average_tau = average_tau;
disc.discrete.rel_amplitude = rel_amplitude;
disc.discrete.abs_max_amplitude = abs_max_amplitude;
disc.discrete.neglected = FOkin.Neglected;
disc.discrete.expfit = FOkin.Expfit;
end

%%
function [mask_lo, mask_common, mask_hi] = overlap(y1, y2)
arguments
    y1 double
    y2 double
end
    mask_lo = (y1 < min(y2));
    mask_common = ((y1 >= min(y2)) & (y1 <= max(y2)));
    mask_hi = (y1 > max(y2));
    if ~any(mask_hi)
        mask_hi = (y2 > max(y1));
    end
end