from django.shortcuts import HttpResponse,render

def hello(request):
    return HttpResponse("Hello World!")

def cat(request):
    
    return render(request,"cat.html")