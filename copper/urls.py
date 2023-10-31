from django.urls import path
from . import views


urlpatterns = [
    path('', views.copper, name='copper')
]
