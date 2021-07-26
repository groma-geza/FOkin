function show_table(self, type)
arguments
    self
    type {mustBeMember(type, {'discrete', 'expfit'})}
end

disc = self.discrete;
if strcmp(type, 'discrete')
    factor = disc.abs_max_amplitude;
    tau_by_ampl = disc.average_tau;
    ampl_by_ampl = disc.rel_amplitude*factor;
    txt = 'Discretized';
    tau_name = 'average tau (';
else
    expfit = disc.expfit;
    factor = expfit.abs_max_amplitude;
    tau_by_ampl = expfit.tau;
    ampl_by_ampl = expfit.rel_amplitude*factor;
    txt = 'Expfit';
    tau_name = 'tau (';
    
end

[tau_by_tau, idx] = sort(tau_by_ampl, 'ascend');
ampl_by_tau = ampl_by_ampl(idx);
control.factor = factor;
control.absolute = true;

if (self.norm == 1)
    normalized = false;
else
    normalized = true;
end

huif= uifigure('Position', [935   277   300   394]);
if normalized
    txt = [txt, ' (nomalized)'];
end
huif.Name = txt; 
huit = uitable(huif, 'Position',  [0, 0, 250, 390]);
huit.Data = [tau_by_ampl, ampl_by_ampl];
huit.ColumnName = {[tau_name, self.options.time_unit, ')'],...
    'amplitude (abs)'};
huit.ColumnFormat = {'shortE', 'short'};
huit.UserData = control;
style = uistyle;
style.HorizontalAlignment = 'center';
addStyle(huit, style)

hmenu1 = uimenu(huif, 'Text', 'Sort');
hitem11 = uimenu(hmenu1, 'Text', 'by tau');
hitem11.MenuSelectedFcn = @MenuSelectedFcn11;
    function MenuSelectedFcn11(~, ~)
        huit.Data = [tau_by_tau, ampl_by_tau];
    end

hitem12 = uimenu(hmenu1, 'Text', 'by amplitude');
hitem12.MenuSelectedFcn = @MenuSelectedFcn12;
    function MenuSelectedFcn12(~, ~)
        huit.Data = [tau_by_ampl, ampl_by_ampl];
    end

hmenu2 = uimenu(huif, 'Text', 'Amplitude');
hitem21 = uimenu(hmenu2, 'Text', 'absolute');
hitem21.MenuSelectedFcn = @MenuSelectedFcn21;
    function MenuSelectedFcn21(~, ~)
        if ~huit.UserData.absolute
            huit.Data(:,2) = huit.Data(:,2)*huit.UserData.factor;
            huit.UserData.absolute = true;
            huit.ColumnName{2} = 'amplitude (abs)';
        end
    end

hitem22 = uimenu(hmenu2, 'Text', 'relative');
hitem22.MenuSelectedFcn = @MenuSelectedFcn22;
    function MenuSelectedFcn22(~, ~)
        if huit.UserData.absolute
            huit.Data(:,2) = huit.Data(:,2)/huit.UserData.factor;
            huit.UserData.absolute = false;
            huit.ColumnName{2} = 'amplitude (rel)';
        end
    end

hmenu3 = uimenu(huif, 'Text', 'Copy');
hitem31 = uimenu(hmenu3, 'Text', 'column of tau');
hitem31.MenuSelectedFcn = @MenuSelectedFcn31;
    function MenuSelectedFcn31(~, ~)
        ctext = sprintf('%.15g\n', huit.Data(:,1));
        clipboard('copy', ctext(1:(end-1)));
    end

hitem32 = uimenu(hmenu3, 'Text', 'column of amplitude');
hitem32.MenuSelectedFcn = @MenuSelectedFcn32;
    function MenuSelectedFcn32(~, ~)
        ctext = sprintf('%.15g\n', huit.Data(:,2));
        clipboard('copy', ctext(1:(end-1)));
    end

hitem33 = uimenu(hmenu3, 'Text', 'both columns');
hitem33.MenuSelectedFcn = @MenuSelectedFcn33;
    function MenuSelectedFcn33(~, ~)
        ctext = sprintf('%.15g\t%.15g\n', huit.Data');
        clipboard('copy', ctext(1:(end-1)));
    end
end