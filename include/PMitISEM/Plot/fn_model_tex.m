function [model_tex, ML] = fn_model_tex(model)
    switch model
        case 't_gas'
            model_tex = '\\textbf{GAS(1,1)-$t$}';          
            ML = false;
        case 't_gas_ML'
            model_tex = '\\textbf{GAS(1,1)-$t$}';   
            ML = true;
        case 't_garch2_noS'
            model_tex = '\\textbf{GARCH(1,1)-$t$}';
            ML = false; 
        case 't_garch2_noS_ML'
            model_tex = '\\textbf{GARCH(1,1)-$t$}';
            ML = true;             
        case 'arch'
            model_tex = '\\textbf{ARCH(1)}';
            ML = false;     
        case 'arch_ML'
            model_tex = '\\textbf{ARCH(1)}';
            ML = true;     
        case 'WN'
            model_tex = '\\textbf{White Noise}';
            ML = false;        
        case 'WN_ML'
            model_tex = '\\textbf{White Noise}';    
            ML = true;        
    end
end