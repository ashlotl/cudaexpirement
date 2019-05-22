import pygame, socket, sys, json, time, math

pygame.init()
screen=pygame.display.set_mode((500,200),pygame.RESIZABLE)
pygame.display.set_caption("Scary port v2084")
UDP_IP = "localhost"
UDP_PORT = 2084

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.connect((UDP_IP, UDP_PORT))

res=int(sys.argv[1])
alt=float(sys.argv[2])
lh=.5

numSamples=0;
ta=lh
while True:
    if ta>=alt:
        break
    #You may notice we project a cube into a sphere. Primarily I'm too lazy to generate something like an icosahedron, secondarily it would not have fine resolution (only subdivisions) if we used an icosahedron, tertiarily this allows us to check for seams, which are screaming alarms something is wrong, and fourth uv spheres get ridiculously fine detail only at the poles, which is a waste.
    numSamples+=res**3-(res-2)**3;
    ta+=lh

points = [[0.0 for i in range(3)] for j in range(numSamples)]
print(len(points))
ta=lh
iterator=0
while (True):
    if ta>=alt:
        break
    xc=-1
    while (True):
        if xc>1:
            break
        yc=-1
        while (True):
            if yc>1:
                break
            zc=-1
            while (True):
                if zc>1:
                    break
                print(xc,yc,zc,ta)
                if xc==-1 or xc==1 or yc==-1 or yc==1 or zc==-1 or zc==1:
                    altAdjust=ta/(xc**2+yc**2.0+zc**2.0)**(1/2)
                    print(ta,xc,yc,zc,iterator,alt)
                    points[iterator]=[xc*altAdjust,yc*altAdjust,zc*altAdjust]
                    iterator+=1
                zc+=2/(res-1)
            yc+=2/(res-1)
        xc+=2/(res-1)
    ta+=lh

sock.send("Arbitrary\n")
largest=0.0
scale=10
speed=1.1
while True:
    changed = False
    for event in pygame.event.get():
        if event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 4:
                print("zooooom...")
                scale/=speed
            if event.button == 5: scale*=speed
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit(0)
        elif event.type == pygame.VIDEORESIZE:
            scrsize = event.size
            width   = event.w
            hight   = event.h
            screen = pygame.display.set_mode(scrsize,pygame.RESIZABLE)
            changed = True
    # print("Receiving...")
    # time.sleep(1)
    print("Waiting...")
    data, addr = sock.recvfrom(numSamples*32)
    data=data.rstrip('\x00')
    print("Parsing...",data)
    if data!="1" and data!="2":
        data=data[0:data.index(']')+1]
        print("new data",data)
        screen.fill((0,0,0))
        pointdata=(json.loads(data))
        itr=0
        w, h = pygame.display.get_surface().get_size()
        for dp in pointdata:
            if dp>largest:
                largest=dp
            # print(int(100*dp/largest))
            color=(float(255*dp/largest),float(100),float(100),float(1))
            pygame.draw.circle(screen,color,(int((points[itr][0]-points[itr][2]+alt)/2*scale+w/2),int((points[itr][1]-points[itr][2]+alt)/2*scale+h/2)),2,0)
            itr+=1
        pygame.display.flip()
