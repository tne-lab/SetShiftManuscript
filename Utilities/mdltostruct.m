function mdlS = mdltostruct(mdl)
    mdlS = struct;
    mdlS.CoefficientNames = mdl.CoefficientNames;
    mdlS.Coefficients = dataset2table(mdl.Coefficients);
    mdlS.Coefficients.Properties.RowNames = mdlS.Coefficients.Name;
    mdlS.NumCoefficients = mdl.NumCoefficients;
    mdlS.ResponseName = mdl.ResponseName;
end