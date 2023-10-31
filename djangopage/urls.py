"""djangopage URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.conf.urls.static import static
from django.conf.urls import url
from django.conf import settings
from django.contrib import admin
from django.urls import path, include
from . import nickelammonia, nickelpure, nickelcitrate, copper, nickelca


urlpatterns = [
    path('admin/', admin.site.urls),
    path('nickelammonia/', include('nickelammonia.urls')),
    path('nickelcitrate/', include('nickelcitrate.urls')),
    path('nickelpure/', include('nickelpure.urls')),
    path('copper/', include('copper.urls')),
    path('nickelca/', include('nickelca.urls')),
    path('', include('home.urls')),
    path('django_plotly_dash/', include('django_plotly_dash.urls')),
    # url(r'^\.well-known/', include('letsencrypt.urls')),
]

if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
