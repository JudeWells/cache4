def get_label(d):
    if d=="101-1000":
        return(550.5)
    if d=="1001-5000":
        return(3000.5)
    if d==">5000":
        return(5000)
    if d=="<100":
        return(100)
    if d=="101-300":
        return(200.5)
    if d=="301-1000":
        return(650.5)
    if d=="1001-3000":
        return(2000.5)
    if d=="3001-10000":
        return(6500.5)
    if d==">10000":
        return(10000) 