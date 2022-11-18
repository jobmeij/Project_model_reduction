function T = CalculateOutput(model,phi,sys)
T = phi.xy*sys.a;
T = T + model.Tamb;
end