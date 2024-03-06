function DbgWaitBar(debug, txt, percentage)
persistent hwb;

if debug == 1
    if percentage == 0
        hwb = waitbar(0, txt);
    elseif percentage == 1
        delete(hwb);
    else
        waitbar(percentage, hwb);
    end
end
end
