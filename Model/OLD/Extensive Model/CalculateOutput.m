function T = CalculateOutput(model,phi,sys)
T = phi.xynorm*sys.a;
T = T + model.Tamb;
end